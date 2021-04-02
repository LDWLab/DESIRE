import io, json, os
from django.http import JsonResponse, HttpResponse, HttpResponseServerError
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.mmcifio import MMCIFIO

def handleCustomUploadStructure (request, strucID):
    '''We will POST all structure chains we need with uniqueIDs.
    Then when we GET them we can list the strucIDs separated by coma 
    this would mean "combine these in one CIF and return them".
    '''
    strucID = strucID.upper()
    if request.method == 'POST':
        from urllib.request import urlopen
        try:
            entities = request.POST.get("entities")
        except:
            return HttpResponseServerError("POST was sent without entities to parse!")
        deStrEnt = json.loads(entities)
        if strucID == "CUST":
            strucString = parseCustomPDB(deStrEnt["stringData"])
            topology = handleTopologyBuilding(deStrEnt["stringData"], "/home/Desire-Server/proorigami-cde-package/cde-root/home/proorigami/")
            request.session[f'TOPOLOGY-{strucID}-{deStrEnt["entityID"]}-{deStrEnt["chainID"]}'] = topology
            #make topology data
            request.session[f'{strucID}-{deStrEnt["entityID"]}-{deStrEnt["chainID"]}'] = strucString
            return JsonResponse("Success!", safe=False)
        for entry in deStrEnt:
            if request.session.get(f'{strucID}-{entry["entityID"]}-{entry["chainID"]}'):
                continue
            ebiURL = f'https://coords.litemol.org/{strucID.lower()}/chains?entityId={entry["entityID"]}&authAsymId={entry["chainID"]}&atomSitesOnly=1'
            #ebiURL = f'https://www.ebi.ac.uk/pdbe/coordinates/{strucID.lower()}/chains?entityId={entityId}&atomSitesOnly=1'
            try:
                data = urlopen(ebiURL)
            except:
                return HttpResponseServerError(f'Failed to fetch coordinates from litemol for PDB {strucID} and entityID {entry["entityID"]} and chain id {entry["chainID"]}.')
            try:
                tempStrucStr = str()
                for line in data:
                    tempStrucStr+=line.decode('UTF-8')
            except:
                return HttpResponseServerError("Failed to parse the provided structure!")
            request.session[f'{strucID}-{entry["entityID"]}-{entry["chainID"]}'] = tempStrucStr
        return JsonResponse("Success!", safe=False)
    
    if request.method == 'GET':
        if request.session.get(strucID):
            stringStruc = request.session[strucID]
            return HttpResponse(stringStruc, content_type="text/plain")
        structureList = list()
        chainToEntity = dict()
        for singleID in strucID.split(','):
            pdbAndEntityAndChain = singleID.split('-')
            try:
                stringData = request.session[singleID]
            except:
                return HttpResponseServerError(f"Requested structure {singleID} was not present in the session! Wait for POST to finish.")
            strucObj = parse_string_structure(stringData, pdbAndEntityAndChain[0])
            chainToEntity[pdbAndEntityAndChain[2]] = pdbAndEntityAndChain[1]
            structureList.append(strucObj)
        try:
            if len(structureList) > 1:
                structureList = combineChainsInSingleStruc(structureList)
            stringStruc = strucToString(structureList[0])
            listStruc = stringStruc.split('\n')
            tempStruc = listStruc[:21]
            for row in listStruc[21:]:
                rowList = row.split()
                if len(rowList) > 10:
                    rowList[7] = chainToEntity[rowList[16]]
                tempStruc.append(' '.join(rowList))
            request.session[strucID] = '\n'.join(tempStruc)
            stringStruc = '\n'.join(tempStruc)
        except:
            return HttpResponseServerError("Failed to parse structures!")
        return HttpResponse(stringStruc, content_type="text/plain")

def parse_string_structure(stringData, strucID):

    parser = MMCIFParser()
    strucFile = io.StringIO(stringData)
    structureObj = parser.get_structure(strucID,strucFile)
    return structureObj

def strucToString(strucObj):
    
    strucFile = io.StringIO("")
    mmCIFio=MMCIFIO()
    mmCIFio.set_structure(strucObj)
    mmCIFio.save(strucFile)
    return strucFile.getvalue()

def combineChainsInSingleStruc(structureList):
    for singleStruc in structureList:
        chains = list(singleStruc.get_chains())
        for removeChain in chains[1:]:
            singleStruc[0].detach_child(removeChain.id)
    for strucToMerge in structureList[1:]:
        chain = list(strucToMerge.get_chains())
        structureList[0][0].add(chain[0])
    return structureList

def parseCustomPDB(stringData):
    parser = PDBParser()
    strucFile = io.StringIO(stringData)
    structureObj = parser.get_structure("CUST",strucFile)
    return strucToString(structureObj)

def handleTopologyBuilding(pdbString, proorigamiLocation):
    from subprocess import Popen, PIPE
    from os import remove, path
    import datetime

    cwd = os.getcwd()
    now = datetime.datetime.now()
    
    fileNameSuffix = "_" + str(now.year) + "_" + str(now.month) + "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) + "_" + str(now.second) + "_" + str(now.microsecond)
    fileLoc = f"{proorigamiLocation}CUSTOMPDB{fileNameSuffix}"
    tempfiles = [f"{fileLoc}.pdb", f"{fileLoc}.svg", f"{fileLoc}.png"]
    for tempf in tempfiles:
        if path.isfile(tempf):
            remove(tempf)
    
    fh = open(f"{fileLoc}.pdb", "w")
    fh.write(pdbString)
    fh.close()

    os.chdir(proorigamiLocation)
    pipe = Popen(f"./make_cartoon.sh.cde {fileLoc}.pdb ; cat {fileLoc}.svg", stdout=PIPE, shell=True)
    output = pipe.communicate()[0]
    os.chdir(cwd)

    if len(output.decode("ascii")) <= 0:
        for removeFile in tempfiles:
            remove(removeFile)
        return HttpResponseServerError("Failed creating topology diagram!\nTry a different structure.")

    svgData = output.decode("ascii")
    for removeFile in tempfiles:
        remove(removeFile)

    return svgData

def getTopology (request, topID):
    if request.session.get(topID):
        topology = request.session[topID]
        return HttpResponse(topology, content_type="text/plain")