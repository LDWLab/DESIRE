import {loadAlignmentViewer} from './loadAlignmentViewer.js'
export function getStructMappingAndTWC (fasta, struc_id, startIndex, stopIndex, ebi_sequence, vueObj){
    if (vm.fasta_data){
        let cleanFasta = vm.fasta_data.replace(/^>Structure sequence\n(.+\n)+?>/i, ">");
        vm.fasta_data = cleanFasta;
    }
    ajax('/mapSeqAln/', {fasta, struc_id}).then(structMappingAndData=>{
        var struct_mapping = structMappingAndData["structureMapping"];
        var largestKey = Math.max(...Object.values(struct_mapping).filter(a=>typeof(a)=="number"))
        var smallestKey = Math.min(...Object.values(struct_mapping).filter(a=>typeof(a)=="number"))
        if ((largestKey != stopIndex || smallestKey != startIndex)&&ebi_sequence){
            ajax('/mapSeqAlnOrig/', {fasta, ebi_sequence, startIndex:1}).then(origStructMappingAndData=>{
                var orig_struct_mapping = origStructMappingAndData["structureMapping"];
                if (structMappingAndData["gapsInStruc"]&&structMappingAndData["gapsInStruc"].length > 0){
                    structMappingAndData["gapsInStruc"].forEach(function(gapTup){
                        let lowMiss = Number(_.invert(struct_mapping)[gapTup[0]]);
                        let topMiss = Number(_.invert(struct_mapping)[gapTup[1]]);
                        if (topMiss-lowMiss>1){
                            for (let i = lowMiss+1; i < topMiss; i++){
                                delete orig_struct_mapping[i];
                            }
                        }
                    })
                }
                assignColorsAndStrucMappings(vueObj, origStructMappingAndData);
            })
        } else {
            assignColorsAndStrucMappings(vueObj, structMappingAndData)
        }
    }).catch(error => {
        var topview = document.querySelector('#topview');
        console.log(error);
        vueObj.topology_loaded = 'error';
        topview.innerHTML = "Failed to load the alignment-structure mapping!<br>Try another structure."
    });
}

var assignColorsAndStrucMappings = function (vueObj, struct_mapping){
    vueObj.poor_structure_map = struct_mapping['BadMappingPositions'];
    vueObj.fasta_data = struct_mapping["amendedAln"];
    vueObj.structure_mapping = struct_mapping["structureMapping"];
    loadAlignmentViewer (vueObj.fasta_data);
    var mapped_aa_properties = mapAAProps(vueObj.aa_properties, vueObj.structure_mapping);
    if (((vueObj.tax_id != null && vueObj.tax_id.length == 2) || (vueObj.custom_aln_twc_flag != null && vueObj.custom_aln_twc_flag == true) || (vueObj.type_tree == 'para'))) {
        if (vueObj.unmappedTWCdata) {
            mapTWCdata(vueObj.structure_mapping, vueObj.unmappedTWCdata, mapped_aa_properties);
        }
    }
    window.mapped_aa_properties = mapped_aa_properties;
    retry(delayedMapping, 10, 1000);
}

var delayedMapping = function (){
    viewerInstanceTop.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
    viewerInstanceTop.pluginInstance.createDomainDropdown();
}

function retry (fn, maxAttempts = 1, delay = 0, attempts = 0) {
    return Promise.resolve()
      .then(fn)
      .catch(err => {
        if (attempts < maxAttempts) {
          return retry (fn, maxAttempts, delay, attempts + 1)
        }
        throw err
      })
  }