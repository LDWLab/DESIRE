import contextlib
import re, os, warnings, io, base64, json
import datetime
import urllib.request
from subprocess import Popen, PIPE
from Bio import AlignIO, BiopythonDeprecationWarning
from io import StringIO

from django.shortcuts import render
from django.http import HttpResponse, JsonResponse, HttpResponseServerError, HttpResponseForbidden
from django.urls import reverse
from django.contrib.sites.shortcuts import get_current_site
import requests

from alignments.models import *
from alignments.taxonomy_views import *
from alignments.residue_api import *
from alignments.structure_api import *
from alignments.fold_api import *
import alignments.alignment_query_and_build as aqab
from TwinCons.bin.TwinCons import slice_by_name
from django.db import connection
from urllib.parse import urlparse

def trim_alignment(concat_fasta, filter_strain):
    '''Reads a fasta string into alignment and trims it down by filter sequence'''
    from alignments.Shannon import species_index_to_aln_index, truncate_aln
    from Bio import AlignIO
    from io import StringIO
    alignment = list(AlignIO.parse(StringIO(concat_fasta), 'fasta'))[0]
    aln_anchor_map, anchor_ix_in_alignment = species_index_to_aln_index(alignment, filter_strain)
    alignment = truncate_aln(alignment, list(aln_anchor_map.keys()), aln_anchor_map=aln_anchor_map)
    return alignment

def calculate_twincons(alignment):
    '''Calculates twincons score given an alignment object.
    Returns data in a list format for the topology viewer'''
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonDeprecationWarning)
    from TwinCons.bin import TwinCons
    list_for_phymeas = ['-as',alignment, '-r', '-mx', 'blosum62']
    alnindex_score, sliced_alns, number_of_aligned_positions, gp_mapping = TwinCons.main(list_for_phymeas)
    list_for_topology_viewer = []
    for alnindex in alnindex_score:
        list_for_topology_viewer.append([alnindex,alnindex_score[alnindex][0]])
    return list_for_topology_viewer

def upload_custom_data_for_mapping(request):
    if request.method == 'POST' and 'filename' in request.FILES:
        data_pairs = []
        file = request.FILES['filename']
        file_iterator = iter(file)
        while True:
            try:
                entry = file_iterator.__next__().decode().strip().split(',')
                data_pairs.append((int(entry[0]), float(entry[1])))
            except StopIteration:
                break
        request.session['csv'] = data_pairs
    if request.method == 'GET':
        data_pairs = request.session.get('csv')
        return JsonResponse(data_pairs, safe = False)

def api_twc_with_upload(request, anchor_structure):
    #### _____________Transform PDBID to taxid______________ ####
    anchor_taxid = pdbid_to_strainid(anchor_structure)

    fastastring = request.session.get('fasta')
    #print('fastastring:\n' + fastastring)

    concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
    #### _____________Trim down the alignment______________ ####
    alignment = trim_alignment(concat_fasta, str(anchor_taxid))

    #### _______________Calculate TwinCons_________________ ####
    list_for_topology_viewer = calculate_twincons(alignment)

    return JsonResponse(list_for_topology_viewer, safe = False)

def constructEbiAlignmentString(fasta, ebi_sequence, startIndex):
    now = datetime.datetime.now()
    fileNameSuffix = "_" + str(now.year) + "_" + str(now.month) + "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) + "_" + str(now.second) + "_" + str(now.microsecond)
    ### BE CAREFUL WHEN MERGING THE FOLLOWING LINES TO PUBLIC; PATHS ARE HARDCODED FOR THE APACHE SERVER ###
    alignmentFileName = "./static/alignment" + fileNameSuffix + ".txt"
    ebiFileName = "./static/ebi_sequence" + fileNameSuffix + ".txt"
    mappingFileName = ebiFileName + ".map"
    fasta = re.sub('>Structure sequence[\s\S]*?>','>',fasta)
    fh = open(alignmentFileName, "w")
    fh.write(fasta)
    fh.close()

    fh = open(ebiFileName, "w")
    fh.write(">Structure sequence\n")
    fh.write(ebi_sequence)
    fh.close()

    shiftIndexBy = 0
    if startIndex > 1:
        shiftIndexBy = startIndex - 1

    pipe = Popen("mafft --quiet --addfull " + ebiFileName + " --mapout " + alignmentFileName + "; cat " + mappingFileName, stdout=PIPE, shell=True)
    output = pipe.communicate()[0]
    decoded_text = output.decode("ascii")
    
    if len(decoded_text) <= 0:
        for removeFile in [alignmentFileName, ebiFileName]:
            os.remove(removeFile)
        return HttpResponseServerError("Failed mapping the polymer sequence to the alignment!\nTry a different structure.")

    text = decoded_text.split('\n#')[1]
    amendedAln = re.sub('>Structure sequence$','',decoded_text.split('\n#')[0])
    outputDict, mapping, firstLine, badMapping = dict(), dict(), True, 0
    for line in text.split('\n'):
        if firstLine:
            firstLine = False
            continue
        row = line.split(', ')
        if len(row) < 3:
            continue
        if row[2] == '-':
            badMapping += 1
            continue
        if row[1] == '-':
            return HttpResponseServerError("Failed mapping the polymer sequence to the alignment!\nTry a different structure.")
        mapping[int(row[2])] = int(row[1]) + shiftIndexBy

    outputDict["structureMapping"] = mapping
    if badMapping > 0:
        outputDict['BadMappingPositions'] = badMapping

    for removeFile in [alignmentFileName, ebiFileName, mappingFileName]:
        os.remove(removeFile)
    outputDict["amendedAln"] = f'>Structure sequence{amendedAln.split(">Structure sequence")[1]}{amendedAln.split(">Structure sequence")[0]}'
    return outputDict

def request_post_data(post_data):
    fasta = post_data["fasta"]
    ebi_sequence = post_data["ebi_sequence"]
    startIndex = int(post_data["startIndex"])
    return fasta, ebi_sequence, startIndex

def make_map_from_alnix_to_sequenceix(request):
    fasta, ebi_sequence, startIndex = request_post_data(request.POST)
    mapping = constructEbiAlignmentString(fasta, ebi_sequence, startIndex)
    if type(mapping) != dict:
        return mapping
    return JsonResponse(mapping, safe = False)

def api_twc_parameterless(request):
    fasta = request.POST["fasta"]
    concat_fasta = re.sub(r'\\n','\n', fasta,flags=re.M)
    list_for_topology_viewer = calculate_twincons(concat_fasta)
    return JsonResponse(list_for_topology_viewer, safe = False)

def api_twc(request, align_name, tax_group1, tax_group2, anchor_structure=''):

    #### _____________Transform PDBID to taxid______________ ####
    if anchor_structure != '':
        anchor_taxid = pdbid_to_strainid(anchor_structure)
        filter_strain = str(Species.objects.filter(strain_id = anchor_taxid)[0].strain).replace(" ", "_")

    fastastring = request.POST.get('fasta')
    if fastastring is None:
        #### _________Query database for the alignment__________ ####
        align_id = Alignment.objects.filter(name = align_name)[0].aln_id

        rawsqls = []
        for parent in [tax_group1, tax_group2]:
            rawsqls.append((aqab.sql_filtered_aln_query(align_id, parent), Taxgroups.objects.get(pk=parent).groupname))
        nogap_tupaln = dict()
        max_alnposition = 0
        for rawsql, parent in rawsqls:
            nogap_tupaln, max_alnposition= aqab.query_to_dict_structure(rawsql, parent, nogap_tupaln, max_alnposition)

        print(nogap_tupaln)
        #### __________________Build alignment__________________ ####
        fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_alnposition)
    
    concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
    #print(concat_fasta)
    
    #### ______________Trim down the alignment______________ ####
    if anchor_structure != '':
        alignment = trim_alignment(concat_fasta, filter_strain)
    else:
        alignment = concat_fasta
    
    #### ________________Calculate TwinCons_________________ ####
    list_for_topology_viewer = calculate_twincons(alignment)
    
    return JsonResponse(list_for_topology_viewer, safe = False)

def minmaxIndex_handler(minIndex, maxIndex):
    if (minIndex == ''):
        minIndex = str(0)
    else:
        minIndex = str(minIndex)
    if (maxIndex == ''):
        maxIndex = str(100000)
    else:
        maxIndex = str(maxIndex)
    return minIndex, maxIndex

def twincons_handler(request, anchor_structure, chain, align_name='', tax_group1='', tax_group2='', minIndex = '', maxIndex = ''):
    from django.urls import resolve
    current_url = resolve(request.path_info).url_name
    minIndex, maxIndex = minmaxIndex_handler(minIndex, maxIndex)
    context = dict()
    context = {
        'pdbid': anchor_structure, 
        'chainid': chain,
        'minIndex' : minIndex,
        'maxIndex' : maxIndex
    }
    if current_url == 'twc_with_upload':
        context['entropy_address'] = "upload/twc-api/"+str(anchor_structure)
    elif current_url == 'custom_csv_data_viewer':
        upload_custom_data_for_mapping(request)
        context['entropy_address'] = "custom-csv-data"
    elif current_url == 'twincons':
        context['entropy_address'] = "twc-api/"+align_name+"/"+str(tax_group1)+"/"+str(tax_group2)+"/"+str(anchor_structure)
    return render(request, 'alignments/twc_detail.html', context)

def entropy(request, align_name, tax_group, anchor_structure):
    from alignments import Shannon
    taxid = pdbid_to_strainid(anchor_structure)
    filter_strain = Species.objects.filter(strain_id = taxid)[0].strain
    align_id = Alignment.objects.filter(name = align_name)[0].aln_id
    polymerid = PolymerData.objects.values("pdata_id").filter(polymeralignments__aln = align_id, strain = taxid)[0]["pdata_id"]
    chainid = Chainlist.objects.values("chainname").filter(polymer = polymerid)[0]["chainname"]
    
    rawsql_result = aqab.sql_filtered_aln_query(align_id, tax_group)
    nogap_tupaln, max_aln_length = aqab.query_to_dict_structure(rawsql_result, Taxgroups.objects.get(pk=tax_group).groupname)
    fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_aln_length)
    #fastastring,max_aln_length = aqab.sql_filtered_aln_query(align_id,tax_group)
    aln_shannon_list = Shannon.main(['-a',fastastring,'-f','fastastring','--return_within','-s',filter_strain])
    #print(aln_shannon_list)
    context = {
        'pdbid': anchor_structure, 
        'chainid': chainid, 
        'shannon_dictionary': aln_shannon_list, 
        'entropy_address':"entropy-api/"+align_name+"/"+str(tax_group)+"/"+str(anchor_structure)
    }
    return render(request, 'alignments/entropy_detail.html', context)

def api_entropy(request, align_name, tax_group, anchor_structure):
    from alignments import Shannon
    import os
    from django.conf import settings
    taxid = pdbid_to_strainid(anchor_structure)
    filter_strain = Species.objects.filter(strain_id = taxid)[0].strain
    align_id = Alignment.objects.filter(name = align_name)[0].aln_id
    
    rawsql_result = aqab.sql_filtered_aln_query(align_id, tax_group)
    nogap_tupaln, max_aln_length = aqab.query_to_dict_structure(rawsql_result, Taxgroups.objects.get(pk=tax_group).groupname)
    fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_aln_length)
    #fastastring,max_aln_length = aqab.sql_filtered_aln_query(align_id,tax_group)
    aln_shannon_list = Shannon.main(['-a',fastastring,'-f','fastastring','--return_within','-s',filter_strain])
    return JsonResponse(aln_shannon_list, safe = False)

def index_test(request):
    some_Alignments = Alignment.objects.all()
    superKingdoms = Taxgroups.objects.raw('SELECT * FROM TaxGroups WHERE\
         TaxGroups.groupLevel = "superkingdom";')
    
    context = {
        'props': list(Taxgroups.objects.values('taxgroup_id', 'groupname')),
        'some_Alignments': some_Alignments,
        'superKingdoms': superKingdoms
    }
    return render(request, 'alignments/index_test.html', context)

def desireAPI (request):
    return render(request, 'alignments/desireAPIindex.html')

def proteinTypes(request):
    if request.method == 'POST' and 'taxIDs' in request.POST:
        taxIDs = request.POST['taxIDs']
        return proteinTypesDirect(request, taxIDs)
    else:
        return []

def allProteinTypes(request):
    allProteinTypes = []
    with connection.cursor() as cursor:
        sql = "select distinct(MoleculeGroup) from Nomenclature order by MoleculeGroup ASC;"
        cursor.execute(sql)
        for row in cursor.fetchall():
            allProteinTypes.append(row[0])
    context = {
        "allProteinTypes" : allProteinTypes
    }
    return JsonResponse(context)

def allSpecies(request):
    allSpecies = []
    with connection.cursor() as cursor:
        sql = "select strain, strain_id from Species order by strain asc;"
        cursor.execute(sql)
        for row in cursor.fetchall():
            allSpecies.append([row[0], row[1]])
    context = {
        "allSpecies" : allSpecies
    }
    return JsonResponse(context)

def proteinTypesDirect(request, concatenatedTaxIds):
    results = []
    with connection.cursor() as cursor:
        if concatenatedTaxIds.endswith(','):
            concatenatedTaxIds = concatenatedTaxIds[:-1]
        for taxID in concatenatedTaxIds.split(','):
            taxID = int(taxID)
            proteinTypesList = []
            # sql = 'SET sql_mode=(SELECT REPLACE(@@sql_mode,\'ONLY_FULL_GROUP_BY\',\'\'));'
            sql = 'select Nomenclature.MoleculeGroup from TaxGroups join Species_TaxGroup on Species_TaxGroup.taxgroup_id = TaxGroups.taxgroup_id join Species on Species.strain_id = Species_TaxGroup.strain_id join Species_Polymer on Species.strain_id = Species_Polymer.strain_id join Polymer_Data on Polymer_Data.strain_id = Species_Polymer.strain_id and Polymer_Data.GI = Species_Polymer.GI and Polymer_Data.nomgd_id = Species_Polymer.nomgd_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id where TaxGroups.taxgroup_id = ' + str(taxID) + ' group by Nomenclature.MoleculeGroup;'

            cursor.execute(sql)
            # results = Taxgroups.objects.raw(sql)
            for row in cursor.fetchall():
                proteinTypesList.append(row[0])
            results.append(proteinTypesList)
    context = {
        'results' : results
    }
    return JsonResponse(context)

def getAlignmentsFilterByProteinTypeAndTaxIds(request):
    if request.method == 'POST' and 'selectedProteinType' in request.POST and 'taxIDs' in request.POST:
        return getAlignmentsFilterByProteinTypeAndTaxIdsDirect(request, request.POST['selectedProteinType'], request.POST['taxIDs'])
    else:
        return []

def getAlignmentsFilterByProteinTypeAndTaxIdsDirect(request, concatenatedProteinTypes, concatenatedTaxIds):
    results = []
    with connection.cursor() as cursor:
        if concatenatedProteinTypes.endswith(','):
            concatenatedProteinTypes = concatenatedProteinTypes[:-1]
        if concatenatedTaxIds.endswith(','):
            concatenatedTaxIds = concatenatedTaxIds[:-1]
        proteinTypes = concatenatedProteinTypes.split(',')
        concatenatedProteinTypes = '\'' + proteinTypes[0] + '\''
        for i in range(1, len(proteinTypes)):
            concatenatedProteinTypes += ', \'' + proteinTypes[i] + '\''
        for taxID in concatenatedTaxIds.split(','):
            alignmentNamesAndPrimaryKeys = []
            sql = "select Alignment.Name, Alignment.Aln_id from Nomenclature join Polymer_Data on Polymer_Data.nomgd_id = Nomenclature.nom_id join Polymer_Alignments on Polymer_Alignments.PData_id = Polymer_Data.PData_id join Alignment on Alignment.Aln_id = Polymer_Alignments.Aln_id join Species_Polymer on Species_Polymer.strain_id = Polymer_Data.strain_id and Species_Polymer.GI = Polymer_Data.GI and Species_Polymer.nomgd_id = Polymer_Data.nomgd_id join Species on Species.strain_id = Species_Polymer.strain_id join Species_TaxGroup on Species.strain_id = Species_TaxGroup.strain_id join TaxGroups on TaxGroups.taxgroup_id = Species_TaxGroup.taxgroup_id where Nomenclature.MoleculeGroup in (" + concatenatedProteinTypes + ") and TaxGroups.taxgroup_id = " + taxID + " group by Alignment.Aln_id;"

            cursor.execute(sql)
            for row in cursor.fetchall():
                alignmentNamesAndPrimaryKeys.append([row[0], row[1]])
            results.append(alignmentNamesAndPrimaryKeys)
    context = {'results' : results}
    return JsonResponse(context)
    # alignmentNamesAndPrimaryKeys = []
    # condition = request.method == 'POST' and 'selectedProteinType' in request.POST and 'taxIDs' in request.POST
    # if condition:
    #     selectedProteinType = request.POST['selectedProteinType']
    #     taxIDs = request.POST['taxIDs']
    #     taxIDsList = []
    #     for taxID in taxIDs:
    #         taxIDsList.append(taxID)
    #     sql = "select Alignment.Name, Alignment.Aln_id from Nomenclature join Polymer_Data on Polymer_Data.nomgd_id = Nomenclature.nom_id join Polymer_Alignments on Polymer_Alignments.PData_id = Polymer_Data.PData_id join Alignment on Alignment.Aln_id = Polymer_Alignments.Aln_id join Species on Species.strain_id = Polymer_Data.strain_id join Species_TaxGroup on Species.strain_id = Species_TaxGroup.strain_id join TaxGroups on TaxGroups.taxgroup_id = Species_TaxGroup.taxgroup_id where Nomenclature.MoleculeGroup = '" + selectedProteinType + "' and TaxGroups.taxgroup_id in (" + str(taxIDsList)[1:-1] + ") group by Alignment.Aln_id;"
    #     print('\n')
    #     print(sql)
    #     print('\n')
    #     # sql = "select Alignment.Name, Alignment.Aln_id from Nomenclature join Polymer_Data on Polymer_Data.nomgd_id = Nomenclature.nom_id join Polymer_Alignments on Polymer_Alignments.PData_id = Polymer_Data.PData_id join Alignment on Alignment.Aln_id = Polymer_Alignments.Aln_id join Species on Species.strain_id = PolymerData.strain_id join Species_TaxGroup on Species_TaxGroup.strain_id = TaxGroups.taxgroup_id where Nomenclature.MoleculeGroup = '" + selectedProteinType + "' group by Alignment.Aln_id and TaxGroups.taxgroup_id in (" + str(taxIDsList)[1:-1] + ");"
    #     with connection.cursor() as cursor:
    #         cursor.execute(sql)
    #         for row in cursor.fetchall():
    #             alignmentNamesAndPrimaryKeys.append([row[0], row[1]])
    # context = {'alignmentNamesAndPrimaryKeys' : alignmentNamesAndPrimaryKeys}
    # return JsonResponse(context)

def getAlignmentsFilterByProteinTypeDirect(request, concatenatedProteinTypes):
    proteinTypes = concatenatedProteinTypes.split(',')
    concatenatedProteinTypes = '\'' + proteinTypes[0] + '\''
    for i in range(1, len(proteinTypes)):
        concatenatedProteinTypes += ', \'' + proteinTypes[i] + '\''
    sql = "select distinct(Name) from Alignment join Polymer_Alignments on Polymer_Alignments.Aln_id = Alignment.Aln_id join Polymer_Data on Polymer_Data.PData_id = Polymer_Alignments.PData_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id where Nomenclature.MoleculeGroup in (" + concatenatedProteinTypes + ");"
    with connection.cursor() as cursor:
        cursor.execute(sql)
        results = []
        for row in cursor.fetchall():
            results.append(row[0])
    context = {
        'results' : results
    }
    return JsonResponse(context)

def getProteinInformationFilterByStrainIDAndProteinNameDirect(request, strain_id, protein_name, internal = False):
    sql = "select Polymer_Data.GI, Polymer_metadata.encoding_location, Polymer_metadata.classification from Alignment join Polymer_Alignments on Alignment.Aln_id = Polymer_Alignments.Aln_id join Polymer_Data on Polymer_Data.PData_id = Polymer_Alignments.PData_id join Species_Polymer on Species_Polymer.nomgd_id = Polymer_Data.nomgd_id and Species_Polymer.GI = Polymer_Data.GI join Species on Species_Polymer.strain_id = Species.strain_id join Polymer_metadata on Polymer_metadata.polymer_id = Polymer_Data.PData_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id where Nomenclature.new_name = '" + protein_name + "' and Species.strain_id = " + str(strain_id)
    with connection.cursor() as cursor:
        cursor.execute(sql)
        results = []
        for row in cursor.fetchall():
            rowDictionary = {
                'GI' : row[0],
                'encoding location' : row[1],
                'classification' : row[2]
            }
            results.append(rowDictionary)
    if internal:
        return results
    context = {
        'results' : results
    }
    return JsonResponse(context)

def getPairwiseAlignmentDirect(request, moleculeType, alignmentName, strainId0, strainId1, internal = False):
    strainIds = [strainId0, strainId1]
    concatenatedStrainIds = '\'' + strainIds[0] + '\''
    for i in range(1, len(strainIds)):
        concatenatedStrainIds += ', \'' + strainIds[i] + '\''
    sql = "select TaxGroups.taxgroup_id, Species.strain from Species join Species_TaxGroup on Species.strain_id = Species_TaxGroup.strain_id join TaxGroups on TaxGroups.taxgroup_id = Species_TaxGroup.taxgroup_id where Species.strain_id in (" + concatenatedStrainIds + ") and TaxGroups.groupLevel = 'genus';"
    genusList = []
    speciesNames = []
    with connection.cursor() as cursor:
        cursor.execute(sql)
        for row in cursor.fetchall():
            genusList.append(row[0])
            speciesNames.append(row[1])
    if len(speciesNames) == 1:
        speciesName0 = speciesName1 = speciesNames[0]
    else:
        speciesName0, speciesName1 = speciesNames
    alignment = string_fasta(request, moleculeType, alignmentName, strainId0 + ',' + strainId1, internal=True)
    alignment = alignment.replace('\\n', '\n')
    # alignment = alignment.replace('\n', '\r')
    pairwiseAlignment = ''
    if (alignment.endswith('\n')) :
        alignment = alignment[0:-1]
        alignmentLines = alignment.splitlines()
        modifiedStrainName0 = speciesName0.replace(' ', '_')
        modifiedStrainName1 = speciesName1.replace(' ', '_')
        for i in range(1, len(alignmentLines), 2):
            titleLine = alignmentLines[i - 1]
            alignmentLine = alignmentLines[i]
            if modifiedStrainName0 in titleLine or modifiedStrainName1 in titleLine:
                pairwiseAlignment += alignmentLine + '\n'
    # return pairwiseAlignment
    if internal:
        return pairwiseAlignment
    return JsonResponse(pairwiseAlignment, safe = False)

def showPairwiseAlignmentDirect(request, moleculeType, alignmentName, strainId0, strainId1):
    pairwiseAlignment = getPairwiseAlignmentDirect(request, moleculeType, alignmentName, strainId0, strainId1, internal = True)
    context = {
        'multipleLinesData' : pairwiseAlignment
    }
    return render(request, 'alignments/multipleLinesVisualizer.html', context)

def showProteinResultsDirect(request, strain_id, gi):
    polymerInformation = getProteinInformationFilterByStrainIDAndProteinNameDirect(request, strain_id, gi, True)
    multipleLinesData = None
    if len(polymerInformation) == 0:
        multipleLinesData = "No results!"
    else:
        polymerInformation = polymerInformation[0]
        multipleLinesData = '\n'.join(['GI: ' + str(polymerInformation['GI']), 'encoding location: ' + str(polymerInformation['encoding location']), 'classification: ' + str(polymerInformation['classification'])])
    context = {
        'multipleLinesData' : multipleLinesData
    }
    return render(request, 'alignments/multipleLinesVisualizer.html', context)

def showStrainInformationFilterByProteinNameDirect(request, protein_name):
    strainInformation = getStrainInformationFilterByProteinNameDirect(request, protein_name, internal = True)
    context = {
        'multipleLinesData' : '\n'.join(str(strainDatum[0]) + ', ' + strainDatum[1] for strainDatum in strainInformation)
    }
    return render(request, 'alignments/multipleLinesVisualizer.html', context)

def showProteinNamesPerStrainIDAndProteinTypeDirect(request, strain_id, proteinType):
    alignmentInformation = getProteinNamesFilterByStrainIDAndProteinTypeDirect(request, strain_id, proteinType, internal = True)
    multipleLinesData = None
    if len(alignmentInformation) == 0:
        multipleLinesData = "No results!"
    else:
        multipleLinesData = '\n'.join(list(foo[0] + ', ' + foo[1] for foo in alignmentInformation))
    context = {
        'multipleLinesData' : multipleLinesData
    }
    return render(request, 'alignments/multipleLinesVisualizer.html', context)

def showAlignmentDirect(request, protein_type, aln_name, tax_group):
    alignment = string_fasta(request, protein_type, aln_name, tax_group, True)
    context = {
        'multipleLinesData' : alignment.replace('\\n', '\n')
    }
    return render(request, 'alignments/multipleLinesVisualizer.html', context)

def getProteinNamesFilterByStrainIDAndProteinTypeDirect(request, strain_id, proteinType, internal = False):
    sql = "select Nomenclature.new_name, Polymer_metadata.Fullseq from Polymer_Data join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id join Species_Polymer on Species_Polymer.GI = Polymer_Data.GI and Species_Polymer.nomgd_id = Polymer_Data.nomgd_id join Species on Species.strain_id = Species_Polymer.strain_id join Polymer_Alignments on Polymer_Alignments.PData_id = Polymer_Data.PData_id join Alignment on Alignment.Aln_id = Polymer_Alignments.Aln_id join Polymer_metadata on Polymer_metadata.polymer_id = Polymer_Data.PData_id where Species.strain_id = " + str(strain_id) + " and Nomenclature.MoleculeGroup = '" + proteinType + "' group by Nomenclature.new_name, Polymer_metadata.Fullseq order by Nomenclature.new_name asc"
    results = []
    with connection.cursor() as cursor:
        cursor.execute(sql)
        for row in cursor.fetchall():
            results.append([row[0], row[1]])
    if internal:
        return results
    context = {
        'results' : results
    }
    return JsonResponse(context)

def getProteinNamesFilterByProteinTypeDirect(request, proteinType):
    sql = "select Nomenclature.new_name from Nomenclature where Nomenclature.MoleculeGroup in ('" + proteinType + "') order by Nomenclature.new_name"
    # sql = "select Polymer_Data.GI from Polymer_Data join Nomenclature on Polymer_Data.nomgd_id = Nomenclature.nom_id where Nomenclature.MoleculeGroup = '" + proteinType + "' order by Polymer_Data.GI asc;"
    results = []
    with connection.cursor() as cursor:
        cursor.execute(sql)
        for row in cursor.fetchall():
            results.append(row[0])
    context = {
        'results' : results
    }
    return JsonResponse(context)

def getStrainInformationFilterByProteinNameDirect(request, aln_name, internal = False):
    sql = "select distinct(Species.strain_id), strain from Species join Species_Polymer on Species_Polymer.strain_id = Species.strain_id join Polymer_Data on Species_Polymer.GI = Polymer_Data.GI and Species_Polymer.nomgd_id = Polymer_Data.nomgd_id join Polymer_Alignments on Polymer_Alignments.PData_id = Polymer_Data.PData_id join Alignment on Alignment.Aln_id = Polymer_Alignments.Aln_id where Alignment.Name = '" + aln_name + "'"
    results = []
    with connection.cursor() as cursor:
        cursor.execute(sql)
        for row in cursor.fetchall():
            results.append([row[0], row[1]])
    if internal:
        return results
    context = {
        'results' : results
    }
    return JsonResponse(context)

def getStrainsFilterByMoleculeGroupAndAlignmentDirect(request, concatenatedMoleculeGroups, concatenatedAlignmentNames):
    moleculeGroups = concatenatedMoleculeGroups.split(',')
    concatenatedMoleculeGroups = '\'' + moleculeGroups[0] + '\''
    for i in range(1, len(moleculeGroups)):
        concatenatedMoleculeGroups += ', \'' + moleculeGroups[i] + '\''
    alignmentNames = concatenatedAlignmentNames.split(',')
    concatenatedAlignmentNames = '\'' + alignmentNames[0] + '\''
    for i in range(1, len(alignmentNames)):
        concatenatedAlignmentNames += ', \'' + alignmentNames[i] + '\''
    sql = 'select Species.strain, Species.strain_id from Species join Species_Polymer on Species_Polymer.strain_id = Species.strain_id join Polymer_Data on Polymer_Data.GI = Species_Polymer.GI and Polymer_Data.nomgd_id = Species_Polymer.nomgd_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id join Polymer_Alignments on Polymer_Alignments.PData_id = Polymer_Data.PData_id join Alignment on Alignment.Aln_id = Polymer_Alignments.Aln_id where Nomenclature.MoleculeGroup in (' + concatenatedMoleculeGroups + ') and Alignment.name in (' + concatenatedAlignmentNames + ');'
    with connection.cursor() as cursor:
        cursor.execute(sql)
        results = []
        for row in cursor.fetchall():
            results.append([row[0], row[1]])
    context = {
        'results' : results
    }
    return JsonResponse(context)

def getAlignmentFilterByNameAndMoleculeGroupTrimByStrainIdDirect(request, concatenatedMoleculeGroups, concatenatedAlignmentNames, concatenatedStrainIds):
    moleculeGroups = concatenatedMoleculeGroups.split(',')
    concatenatedMoleculeGroups = '\'' + moleculeGroups[0] + '\''
    for i in range(1, len(moleculeGroups)):
        concatenatedMoleculeGroups += ', \'' + moleculeGroups[i] + '\''
    alignmentNames = concatenatedAlignmentNames.split(',')
    concatenatedAlignmentNames = '\'' + alignmentNames[0] + '\''
    for i in range(1, len(alignmentNames)):
        concatenatedAlignmentNames += ', \'' + alignmentNames[i] + '\''
    sql = 'select Polymer_metadata.Fullseq, Species.strain from Polymer_metadata join Polymer_Data on Polymer_metadata.polymer_id = Polymer_Data.PData_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id join Species_Polymer on Polymer_Data.GI = Species_Polymer.GI and Polymer_Data.nomgd_id = Species_Polymer.nomgd_id join Species on Species.strain_id = Species_Polymer.strain_id join Polymer_Alignments on Polymer_Alignments.PData_id = Polymer_Data.PData_id join Alignment on Polymer_Alignments.Aln_id = Alignment.Aln_id where Nomenclature.MoleculeGroup in (' + concatenatedMoleculeGroups + ') and Alignment.Name in (' + concatenatedAlignmentNames + ') and Species.strain_id in (' + concatenatedStrainIds + ');'
    with connection.cursor() as cursor:
        cursor.execute(sql)
        results = []
        for row in cursor.fetchall():
            results.append(row[0])
    context = {
        'results' : results
    }
    return JsonResponse(context)

def index_orthologs(request):
    if request.method == 'GET' and 'custom_propensity_data' in request.FILES:
        propensity_indices_file = request.FILES['custom_propensity_data']
        propensity_indices_string = ''
        for propensity_part in propensity_indices_file.chunks():
            propensity_indices_string += propensity_part.decode()
    # if request.method == 'POST' and 'custom_propensity_data' in request.FILES:
    #     propensity_indices_file = request.FILES['custom_propensity_data']
    #     propensity_indices_string = ''
    #     for propensity_part in propensity_indices_file.chunks():
    #         propensity_indices_string += propensity_part.decode()
    #     print ("propensity_indices_string: " + propensity_indices_string)
    return render(request, 'alignments/index_orthologs.html')
def index(request):
    return render(request, 'alignments/index.html')

def visualizer(request, align_name, tax_group1, tax_group2, anchor_structure = ''):
    twc_api_url = "http://127.0.0.1:8000/orthologs/twc-api/" + align_name + "/" + str(tax_group1) + "/" + str(tax_group2) + "/" + anchor_structure
    context = {
        "twc_api_url" : twc_api_url
    }
    return render(request, 'alignments/simpleVisualization.html', context)
    #return visualizerHelper(request, align_name + "/" + str(tax_group1) + "/" + str(tax_group2) + "/" + anchor_structure)

def upload_custom_data(request):
    return render(request, 'alignments/upload_custom_data.html')

def paralog_entry_form(request):
    return render(request, 'alignments/index_paralogs.html')

def paralog_display_entropy(request, align_name, fold1, fold2):
    pass
    #return render(request, 'alignments/twc_detail.html', context)

def extract_species_list(fastastring):
    '''Filters out species from a fastastring'''
    unf_species_list = [x.split('\\')[0] for x in fastastring.split('>')[1:]]
    filtered_spec_list = [re.sub('_',' ', re.sub(r'^.*?_', '', x)) for x in unf_species_list]
    return filtered_spec_list

def extract_gap_only_cols(fastastring):
    '''Extracts positions in the fastastring that are only gaps'''
    unf_seq_list = [x.split('\\n')[1] for x in fastastring.split('>')[1:]]
    list_for_intersect = list()
    for sequence in unf_seq_list:
        iterator = re.finditer('-', sequence)
        gap_positions = [m.start(0) for m in iterator]
        list_for_intersect.append(gap_positions)
    gap_only_cols = list(set(list_for_intersect[0]).intersection(*list_for_intersect))
    return gap_only_cols

def construct_dict_for_json_response(response_data):
    '''Takes list of datas for Json response.
    Check their types and assigns names to each.
    Returns them as a dictionary.'''
    response_dict = dict()
    for entry in response_data:
        if type(entry) == str:
            response_dict['Alignment'] = entry
            continue
        if type(entry) == bool:
            response_dict['TwinCons'] = entry
            continue
        if type(entry) == list:
            if all(isinstance(item, int) for item in entry):
                response_dict['Gap-only columns'] = entry
                continue
            if all(isinstance(item, list) for item in entry):
                response_dict['AA frequencies'] = entry
                continue
            if all(isinstance(item, str) for item in entry):
                response_dict['Sequence names'] = entry
                continue
    return response_dict

def simple_fasta(request, aln_id, tax_group, internal=False):
    rawsqls = []
    if type(tax_group) == int:
        tax_group = str(tax_group)
    for parent in tax_group.split(','):
        rawsqls.append((aqab.sql_filtered_aln_query(aln_id, parent), Taxgroups.objects.get(pk=parent).groupname))

    nogap_tupaln = dict()
    max_alnposition = 0

    for rawsql, parent in rawsqls:
        nogap_tupaln, max_alnposition= aqab.query_to_dict_structure(rawsql, parent, nogap_tupaln, max_alnposition)
    
    fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_alnposition)
    
    if internal:
        return fastastring
    
    concat_fasta, twc, gap_only_cols, filtered_spec_list, alignment_obj = calculateFastaProps(fastastring)
    response_dict = construct_dict_for_json_response([concat_fasta,filtered_spec_list,gap_only_cols,frequency_list,twc])

    return JsonResponse(response_dict, safe = False)

def string_fasta(request, protein_type, aln_name, tax_group, internal=False):
    if type(tax_group) == int:
        tax_group = str(tax_group)
    elif type(tax_group) == list:
        tax_group = ','.join(tax_group)
    protein_type = protein_type.split(',')
    for i in range(len(protein_type)):
        protein_type[i] = '\'' + protein_type[i] + '\''
    protein_type=','.join(protein_type)
    with connection.cursor() as cursor:
        sql = 'SET sql_mode=(SELECT REPLACE(@@sql_mode,\'ONLY_FULL_GROUP_BY\',\'\'));'
        cursor.execute(sql)
        sql = "select Alignment.*, Polymer_Data.*, Nomenclature.* from Alignment join Polymer_Alignments on Polymer_Alignments.Aln_id = Alignment.Aln_id join Polymer_Data on Polymer_Data.PData_id = Polymer_Alignments.PData_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id join Species on Polymer_Data.strain_id = Species.strain_id join Species_TaxGroup on Species_TaxGroup.strain_id = Species.strain_id join TaxGroups on Species_TaxGroup.taxgroup_id = Species_TaxGroup.taxgroup_id where Alignment.Name = '" + aln_name + "' and MoleculeGroup in (" + protein_type + ") and TaxGroups.taxgroup_id in (" + tax_group + ") group by Alignment.Name;"
        cursor.execute(sql)
        raw_result = aqab.dictfetchall(cursor)
    return simple_fasta(request, raw_result[0]['Aln_id'], tax_group, internal)
    # return JsonResponse("Test", safe=False)
    # from django.db import connection
    # raw_sql = "select * from Nomenclature where MoleculeGroup = '" + protein_type + "' and new_name = '" + aln_name + "';"
    # with connection.cursor() as cursor:
    #     cursor.execute(raw_sql)
    #     raw_result = aqab.dictfetchall(cursor)
    # if len(raw_result) == 0:
    #     raise Http404("We do not have this combination of arguments in our database.")
    
    # rawsqls = []
    
    # nogap_tupaln = dict()
    # max_alnposition = 0

    # for rawsql, parent in rawsqls:
    #     nogap_tupaln, max_alnposition= aqab.query_to_dict_structure(rawsql, parent, nogap_tupaln, max_alnposition)
    
    # fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_alnposition)
    
    # concat_fasta, twc, gap_only_cols, filtered_spec_list, alignment_obj = calculateFastaProps(fastastring)
    # response_dict = construct_dict_for_json_response([concat_fasta,filtered_spec_list,gap_only_cols,frequency_list,twc])

    # return JsonResponse(response_dict, safe = False)

def getGenusFromStrainIdsDirect(request, concatenatedStrainIds):
    strainIds = concatenatedStrainIds.split(',')
    concatenatedStrainIds = '\'' + strainIds[0] + '\''
    for i in range(1, len(strainIds)):
        concatenatedStrainIds += ', \'' + strainIds[i] + '\''
    sql = "select TaxGroups.taxgroup_id from Species join Species_TaxGroup on Species.strain_id = Species_TaxGroup.strain_id join TaxGroups on TaxGroups.taxgroup_id = Species_TaxGroup.taxgroup_id where Species.strain_id in (" + concatenatedStrainIds + ") and TaxGroups.groupLevel = 'genus';"
    genusList = []
    with connection.cursor() as cursor:
        cursor.execute(sql)
        for row in cursor.fetchall():
            genusList.append(row[0])
    context = {
        'genusList' : genusList
    }
    return JsonResponse(context)

def string_fasta_two_strains(request, protein_type, aln_id, strain_id_0, strain_id_1):
    if type(strain_id_0) == int:
        strain_id_0 = str(strain_id_0)
    if type(strain_id_1) == int:
        strain_id_1 = str(strain_id_1)
    protein_type = protein_type.split(',')
    for i in range(len(protein_type)):
        protein_type[i] = '\'' + protein_type[i] + '\''
    protein_type=','.join(protein_type)
    sql = ""
    with connection.cursor() as cursor:
        cursor.execute(sql)

def calculateFastaProps(fastastring):
    concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
    alignment_obj = AlignIO.read(StringIO(concat_fasta), 'fasta')
    twc = False
    if (len(alignment_obj) < 1000):
        sliced_alns = slice_by_name(alignment_obj)
        if len(sliced_alns.keys()) == 2:
            twc = True
    gap_only_cols = extract_gap_only_cols(fastastring)
    filtered_spec_list = extract_species_list(fastastring)
    return concat_fasta, twc, gap_only_cols, filtered_spec_list, alignment_obj

def rProtein(request, align_name, tax_group):
    #if tax_group == 0 - no filter
    align_id = Alignment.objects.filter(name = align_name)[0].aln_id
    fastastring = simple_fasta(request, align_id, tax_group, internal=True)
    print(fastastring)
    #fastastring,max_aln_length = aqab.sql_filtered_aln_query(align_id,tax_group)
    context = {'fastastring': fastastring, 'aln_name':str(Alignment.objects.filter(aln_id = align_id)[0].name)}
    return render(request, 'alignments/detail.html', context)

def rRNA(request, align_name, tax_group):
    from Bio.SeqUtils import IUPACData
    align_id = Alignment.objects.filter(name = align_name)[0].aln_id
    rawsql_result = aqab.sql_filtered_aln_query(align_id, tax_group)
    nogap_tupaln = dict()
    nogap_tupaln, max_aln_length = aqab.query_to_dict_structure(rawsql_result, Taxgroups.objects.get(pk=tax_group).groupname, nogap_tupaln)
    fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_aln_length, IUPACData.unambiguous_rna_letters)
    #fastastring,max_aln_length = aqab.sql_filtered_aln_query(align_id,tax_group)
    context = {'fastastring': fastastring, 'aln_name':str(Alignment.objects.filter(aln_id = align_id)[0].name)}
    return render(request, 'alignments/rRNA.html', context)

def validate_fasta_string(fastaString):
    malicious_strings = [
        "eval\(unescape",
        "base64_decode\(",
        "substr\(md5\(strrev\(",
        "cwd = @getcwd\(\);",
        "chr\(\(ord\(",
        "gzinflate\(base64_decode\(",
        "php_uname\(\)\" \"] = chr\(ord\(",
        "cwd\[strlen\(\$cwd\)",
        "ini_get\('safe_mode'\);",
        "=\"\x62\"",
        "\"+ r + \"&r=\" + document.referrer;\"",
        "if\(strtoupper\(substr\(PHP_OS, 0, 3\) \) == \"WIN\"\)",
        "window.top.location.href=\"http://",
        "@ini_get\(\"disable_functions\"\)",
        "\$g3='';\$g3.=\$r;\$g3.=\$h;\$g3.=\$y",
    ]
    for mstring in malicious_strings:
        regex = re.compile(mstring)
        if re.search(regex, fastaString):
            return False
    return True

# trims fasta by a list of indices
def trim_fasta_by_index(input_file, indices):
    from Bio import AlignIO
    align = AlignIO.read(input_file, "fasta")
    trimmed_align = align[:,int(indices[0]):int(indices[0])] # initialize align object
    for i in indices.split(','):
        if (int(i)-1 < 0):
            continue
        trimmed_align += align[:,int(i)-1:int(i)]
    return trimmed_align

def propensity_data_custom (request):
    response = propensity_data(request, None, None)
    return response

def propensity_data(request, aln_id, tax_group):
    from io import StringIO
    import alignments.propensities as propensities

    if request.method == 'POST' and 'customFasta' in request.POST and aln_id is None:
        fastastring = request.POST['customFasta']
    else:
        fastastring = simple_fasta(request, aln_id, tax_group, internal=True).replace('\\n', '\n')
    fasta = StringIO(fastastring)

    if request.method == 'POST' and 'indices' in request.POST:
        indices = request.POST['indices']
        trimmed_fasta = trim_fasta_by_index(fasta, indices)
        fasta = StringIO(format(trimmed_fasta, 'fasta'))

    aa = propensities.aa_composition(fasta, reduced = False)
    fasta.seek(0) # reload the fasta object
    red_aa = propensities.aa_composition(fasta, reduced = False)

    data = {
        'aln_id' : aln_id,
        'tax_group' : tax_group,
        'reduced alphabet' : red_aa,
        'amino acid' : aa
    }
    return JsonResponse(data)

def permutation_data(request):
    from io import StringIO
    from Bio import AlignIO
    fasta_variable_name = 'customFasta'
    indices_variable_name = 'indices'
    if not (request.method == 'POST' and fasta_variable_name in request.POST and indices_variable_name in request.POST):
        raise KeyError(fasta_variable_name + " and " + indices_variable_name + " not found in request.POST")
    alignment = AlignIO.read(StringIO(request.POST[fasta_variable_name]), "fasta")
    row_dimension = len(alignment)
    column_dimension = alignment.get_alignment_length()
    # Initialize an empty alignment.
    new_alignment = alignment[:, 0:0]
    for msa_index_pair in re.sub(r"\s*", "", request.POST[indices_variable_name]).split(","):
        dash_index = msa_index_pair.find("-")
        if dash_index >= 0:
            # Both a start index and end index are provided.
            # Convert from a 1-based index system to a 0-based index system.
            start_msa_index = int(msa_index_pair[:dash_index]) - 1
            end_msa_index = int(msa_index_pair[dash_index + 1:])
        else:
            # Only one index provided.
            # Convert from a 1-based index system to a 0-based index system.
            end_msa_index = int(msa_index_pair)
            start_msa_index = end_msa_index - 1
        # Append to the new alignment.
        new_alignment += alignment[:, start_msa_index:end_msa_index]

    # Initialize the label suffixes
    label_suffixes = [""] * row_dimension
    for row_index in range(row_dimension):
        base_column_index = start_msa_index
        while base_column_index < column_dimension and alignment[row_index, base_column_index] == '-':
            base_column_index += 1
        previous_column_index_gap_flag = False
        for column_index in range(base_column_index + 1, end_msa_index):
            column_index_gap_flag = alignment[row_index, column_index] == '-'
            if column_index_gap_flag != previous_column_index_gap_flag:
                if column_index_gap_flag:
                    if column_index == base_column_index + 1:
                        label_suffixes[row_index] += str(base_column_index + 1) + ", "
                    else:
                        label_suffixes[row_index] += str(base_column_index + 1) + "-" + str(column_index + 1) + ", "
                else:
                    base_column_index = column_index
            previous_column_index_gap_flag = column_index_gap_flag
        if alignment[row_index, end_msa_index - 1] != '-':
            if base_column_index == end_msa_index - 1:
                label_suffixes[row_index] += str(base_column_index + 1) + ", "
            else:
                label_suffixes[row_index] += str(base_column_index + 1) + "-" + str(end_msa_index) + ", "

    # Label the fasta-file entries.
    for row_index, label_suffix in enumerate(label_suffixes):
        if len(label_suffix) > 0:
            new_alignment[row_index].id = re.sub('([^_]*)_?$', r'\1_', new_alignment[row_index].id) + label_suffix[:-2]
            new_alignment[row_index].description = ""
    permutation_string = StringIO(format(new_alignment, 'fasta')).getvalue()

    now = datetime.datetime.now()
    fileNameSuffix = "_" + str(now.year) + "_" + str(now.month) + "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) + "_" + str(now.second) + "_" + str(now.microsecond)
    write_to_location = "/home/Desire-DEV/PVDev/static/"
    alignmentFilePath = write_to_location + "permuted_alignment" + fileNameSuffix + ".fasta"

    fh = open(alignmentFilePath, "w")
    fh.write(permutation_string)
    fh.close()

    hhalignOutputFilePath = write_to_location + "test.hhr" + fileNameSuffix
    # return JsonResponse("Output path: " + hhalignOutputFilePath, safe = False)
    output = hhalign(alignmentFilePath, "/home/blastdb/ecodf/ecod_version20220113/20220113/e4gvlA4.hhm", hhalignOutputFilePath)

    os.remove(alignmentFilePath)
    os.remove(hhalignOutputFilePath)

    return output

def hhsearch(input_file_path, input_database_path, output_file_path, max_residues=550000, threshold_percentage=None, add_cons_flag=True, parse_output_flag=False):
    if (max_residues <= 0):
        raise ValueError('The input max_residues value must be greater than zero')
    hhsearch_command = '/usr/local/bin/hh-suite/bin/hhsearch -i ' + input_file_path + ' -d ' + input_database_path + ' -o ' + output_file_path  + ' -maxres ' + str(max_residues)
    if (not threshold_percentage is None):
        if (threshold_percentage < 0 or threshold_percentage > 100):
            raise ValueError('The input threshold_percentage value must be between 0 and 100 inclusively.')
        hhsearch_command += ' -M ' + str(threshold_percentage)
    if add_cons_flag:
        hhsearch_command += ' -add_cons'
    os.system(hhsearch_command)
    if (parse_output_flag):
        return parse_hh_output(output_file_path)

def hhalign(input_query_file_path, input_template_file_path, output_file_path):
    hhalign_command = '/usr/local/bin/hh-suite/bin/hhalign -i ' + input_query_file_path + ' -t ' + input_template_file_path + ' -o ' + output_file_path + ' No Hit Prob E-value P-value Score SS Cols'
    os.system(hhalign_command)
    return parse_hh_output(output_file_path)

def parse_hh_output(hh_output_file_path):
    file_handler = open(hh_output_file_path, 'r')
    lines = file_handler.readlines()
    file_handler.close()
    hit_lines = []
    for index, line in enumerate(lines):
        if re.match('\s*No\s+Hit\s+Prob\s+E-value\s+P-value\s+Score\s+SS\s+Cols\s+Query\s+HMM\s+Template\s+HMM\s*', line):
            for line in lines[index + 1:]:
                if re.match('^\s*$', line):
                    # Stop at the empty line.
                    break
                hit_lines.append(line.strip())
            break
    return parse_hh_output_hit_lines(hit_lines)

def parse_hh_output_hit_lines(hh_output_hit_lines):
    parsed_output = []
    pattern = '^\s*(\d+)\s*e([\d\w]{4})([\d\w]+)(\d)\s+\w+:(\d+)\s*-\s*(\d+)((?:\s*,\s*\w:\d+\s*-\s*\d+\s*)*)\s+\d+\s+([\d\.]+)'
    for line_index, line in enumerate(hh_output_hit_lines):
        regex_results = re.search(pattern, line)
        if regex_results is None:
            raise ValueError(f"The hh output line #{line_index} does not match the expected hit format. Its results could not be parsed.")
        regex_results_groups = regex_results.groups()
        start_variable_name = "start"
        stop_variable_name = "stop"
        _range = [{
            start_variable_name: int(regex_results_groups[4]),
            stop_variable_name: int(regex_results_groups[5])
        }]
        extended_range_string = regex_results_groups[6]
        extended_range_string = re.sub("^[\s,]+", "", extended_range_string)
        extended_range_string = re.sub("[\s,]+$", "", extended_range_string)
        if len(extended_range_string) > 0:
            for range_portion in extended_range_string.split(","):
                range_portion_regex = re.search("\s*(\d+)\s*-\s*(\d+)", range_portion)
                range_portion_regex_groups = range_portion_regex.groups()
                _range.append({
                    start_variable_name: int(range_portion_regex_groups[0]),
                    stop_variable_name: int(range_portion_regex_groups[1])
                })
        parsed_output.append({
            "hit_no": int(regex_results_groups[0]),
            "pdb_id": regex_results_groups[1],
            "chain_id": regex_results_groups[2],
            "pseudo_entity_id": int(regex_results_groups[3]),
            "range": _range,
            "match_percentage": float(regex_results_groups[7])
        })
    return JsonResponse(parsed_output, safe=False)

# def parse_hh_output(hh_output_file_path):
#     file_handler = open(hh_output_file_path, 'r')
#     lines = file_handler.readlines()
#     file_handler.close()
#     pattern = re.compile('^\\s*No\\s+Hit\\s+Prob\\s+E-value\\s+P-value\\s+Score\\s+SS\\s+Cols\\s+Query HMM\\s+Template HMM\\s*$')
#     line_index = 0
#     while (not pattern.match(lines[line_index])):
#         line_index += 1
#     title_line = lines[line_index]
#     index_of_start_of_no = title_line.index('No')
#     index_of_start_of_hit = title_line.index('Hit', index_of_start_of_no + len('No'))
#     index_of_start_of_prob = title_line.index('Prob', index_of_start_of_hit + len('Hit'))
#     index_of_start_of_e_value = title_line.index('E-value', index_of_start_of_prob + len('Prob'))
#     index_of_start_of_p_value = title_line.index('P-value', index_of_start_of_e_value + len('E-value'))
#     index_of_start_of_score = title_line.index('Score', index_of_start_of_p_value + len('P-value'))
#     index_of_start_of_ss = title_line.index('SS', index_of_start_of_score + len('Score'))
#     index_of_start_of_cols = title_line.index('Cols', index_of_start_of_ss + len('SS'))
#     index_of_start_of_query_hmm = title_line.index('Query HMM', index_of_start_of_cols + len('Cols'))
#     index_of_start_of_template_hmm = title_line.index('Template HMM', index_of_start_of_query_hmm + len('Query HMM'))
#     # Skip over the matching title line.
#     line_index += 1
#     start_data_lines_index = line_index
#     pattern = re.compile('^\\s*$')
#     while line_index < len(lines) and (not pattern.match(lines[line_index])):
#         line_index += 1
#     data_lines = lines[start_data_lines_index:line_index]
#     mapped_data = []
#     for data_line_index in range(len(data_lines)):
#         data_line = data_lines[data_line_index]
#         no_substring = data_line[index_of_start_of_no - 1 : index_of_start_of_hit - 1].strip()
#         hit_substring = data_line[index_of_start_of_hit - 1 : index_of_start_of_prob - 1].strip()
#         prob_substring = data_line[index_of_start_of_prob - 1 : index_of_start_of_e_value - 1].strip()
#         e_value_substring = data_line[index_of_start_of_e_value - 1 : index_of_start_of_p_value - 1].strip()
#         p_value_substring = data_line[index_of_start_of_p_value - 1 : index_of_start_of_score - 1].strip()
#         score_substring = data_line[index_of_start_of_score - 1 : index_of_start_of_ss - 1].strip()
#         ss_substring = data_line[index_of_start_of_ss - 1 : index_of_start_of_cols - 1].strip()
#         cols_substring = data_line[index_of_start_of_cols - 1 : index_of_start_of_query_hmm - 1].strip()
#         query_hmm_substring = data_line[index_of_start_of_query_hmm - 1 : index_of_start_of_template_hmm - 1].strip()
#         template_hmm_substring = data_line[index_of_start_of_template_hmm - 1 :].strip()
#         mapped_data += [{
#             'No' : int(no_substring),
#             'Hit' : hit_substring,
#             'Prob' : float(prob_substring),
#             'E-value' : float(e_value_substring),
#             'P-value' : float(p_value_substring),
#             'Score' : float(score_substring),
#             'SS' : float(ss_substring),
#             'Cols' : int(cols_substring),
#             'Query HMM' : query_hmm_substring,
#             'Template HMM' : template_hmm_substring,
#         }]
#     return mapped_data

def proxy(request):
    url = request.GET['url']
    domain = urlparse(url).netloc
    domain_list = ['prodata.swmed.edu']
    if domain not in domain_list:
        return HttpResponseForbidden(
            "This proxy is only for " + ", ".join(domain_list) + "."
        )
    try:
        # http://prodata.swmed.edu/ecod/complete/search?kw=4upz&target=domain&page=4
        if domain == 'prodata.swmed.edu':
            url += '&target=' + request.GET['target'] + '&page=' + request.GET['page']
        proxied_response = requests.get(url)
        if proxied_response.status_code == 200:
            response = HttpResponse(proxied_response.text)
            return response
        else:
            raise Http404
    except:
        raise Http404

def propensities(request, align_name, tax_group):
    aln_id = Alignment.objects.filter(name = align_name)[0].aln_id
    propensity_data = reverse('alignments:propensity_data', kwargs={'aln_id': aln_id, 'tax_group' : tax_group})

    names = []
    if type(tax_group) == int:
        tax_group = str(tax_group)
    for group in tax_group.split(','):
        names.append(Taxgroups.objects.get(pk=group).groupname)

    context = {
        "propensity_data" : propensity_data, 
        "align_name" : align_name,
        "tax_name" : ', '.join(names)
    }
    
    return render(request, 'alignments/propensities.html', context)

def flushSession (request):
    try:
        request.session.flush()
    except:
        return HttpResponseServerError ("Failed to flush the session!")
    return HttpResponse ("Success!")

def ecodPassThroughQuery(request):
    '''Request a password protected URL from our website that returns a JSON object.
    '''
    baseURL = 'http://'+get_current_site(request).domain
    url = baseURL+request.GET['url']
    if ('&format=json' not in url):
        url += '&format=json'
    req = urllib.request.Request(url)
    #username = os.environ['DJANGO_USERNAME']
    #password = os.environ['DJANGO_PASSWORD']
    #credentials = (f'{username}:{password}')
    credentials = ('website:desire_RiboVision3')
    encoded_credentials = base64.b64encode(credentials.encode('ascii'))
    req.add_header('Authorization', 'Basic %s' % encoded_credentials.decode("ascii"))

    response = urllib.request.urlopen(req)
    encoding = response.info().get_content_charset('utf-8')
    data = response.read()
    return JsonResponse(json.loads(data.decode(encoding)), safe=False)


def parse_string_structure(stringData, strucID):
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser()
    strucFile = io.StringIO(stringData)
    structureObj = parser.get_structure(strucID,strucFile)
    return structureObj
