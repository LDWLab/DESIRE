function handleMaskingRanges(mask_range){
    vm.masking_range = mask_range;
    window.masking_range_array = null;
    if (isCorrectMask(mask_range)) {   
        var topviewer = document.getElementById("PdbeTopViewer");
        topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);   
        if(window.custom_prop) {
            topviewer.pluginInstance.getAnnotationFromRibovision(window.custom_prop); 
        }
        window.masked_array = initializeMaskedArray();          
        var selectedIndex = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox').selectedIndex;

        var index = 4;
        while(index < topviewer.pluginInstance.domainTypes.length) {
            colorResidue(index, window.masked_array);
            index++;
        }
        let selectedData = topviewer.pluginInstance.domainTypes[selectedIndex]
        
        if (selectedData.data){
            topviewer.pluginInstance.updateTheme(selectedData.data); 
            window.viewerInstance.visual.select({data: selectSections_RV1.get(selectedData.label), nonSelectedColor: {r:255,g:255,b:255}});
            }
        vm.correct_mask = 'True';
    } else {
        vm.correct_mask = 'False';
    }
}
function handleFilterRange(filter_range) {
    const temp_array = filter_range.split('-');
    if (filter_range.match(/^\d+-\d+/) && Number(temp_array[0]) < Number(temp_array[1])) {
        vm.filter_range = filter_range;
        window.filterRange = temp_array.join(",");
        var topviewer = document.getElementById("PdbeTopViewer");
        var selectedIndex = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox').selectedIndex;
        topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);   
        viewerInstance.visual.update({
            customData: {
                url: `https://www.ebi.ac.uk/pdbe/coordinates/${window.pdblower}/residueRange?entityId=${topviewer.entityId}&range=${filter_range}&encoding=bcif`,
                format: 'cif',
                binary:true },
            assemblyId: '1',
            subscribeEvents: true
        });
        viewerInstance.events.loadComplete.subscribe(() => { 
            let selectedData = topviewer.pluginInstance.domainTypes[selectedIndex];
            let select_sections = selectSections_RV1.get(selectedData.label).slice(Number(temp_array[0]), Number(temp_array[1])+1);
            window.viewerInstance.visual.select({
            data: select_sections,
            nonSelectedColor: {r:255,g:255,b:255}});
            //var selectedDomain = topviewer.pluginInstance.domainTypes[selectedIndex];
            //topviewer.updateTheme(selectedDomain.data);
         });
         topviewer.pluginInstance.initPainting(window.select_sections)
         /*let selectedData = topviewer.pluginInstance.domainTypes[selectedIndex];
         topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);   
         topviewer.pluginInstance.updateTheme(selectedData.data); */
    }else{
        //
    }
}
function isCorrectMask(mask_range){
    window.masking_range_array = null;
    if (mask_range.match(/^(\d+-\d+;)+$/)) {
        var temp_array = mask_range.split(';').join('-').split('-');
        temp_array = temp_array.slice(0, -1)
        var i = 0;
        var isCorrect = true;
        while(i < temp_array.length) {
            if(i % 2 == 0) {
                if(Number(temp_array[i]) > Number(temp_array[i + 1])) {
                    isCorrect = false;
                }
            }
            i = i + 1;
        }
        window.masking_range_array = temp_array;
    }
    return isCorrect;
}
function initializeMaskedArray() {
    var topviewer = document.getElementById("PdbeTopViewer");
    var masked_array = [];
    var j = 0;
    while(j < mapped_aa_properties.get(topviewer.pluginInstance.domainTypes[4].label).length) {
        masked_array[j] = false;
        var i = 0;
        while(i < window.masking_range_array.length && !masked_array[j]) {
            if(j >= window.masking_range_array[i] && j <= window.masking_range_array[i + 1]) {
                masked_array[j] = true;
            }
            i = i+2;
        }
        j = j+1;
    }
    return masked_array;
}
function colorResidue(index, masked_array) {
    var topviewer = document.getElementById("PdbeTopViewer");
    var f = 0;
    while(f < topviewer.pluginInstance.domainTypes[4].data.length) {
        if(!masked_array[f] && topviewer.pluginInstance.domainTypes[index].data[f]) {
            topviewer.pluginInstance.domainTypes[index].data[f].color = "rgb(255,255,255)";
            topviewer.pluginInstance.domainTypes[index].data[f].tooltipMsg = "NaN";                   
            selectSections_RV1.get(topviewer.pluginInstance.domainTypes[index].label)[f].color = {r: 255, g: 255, b: 255};

        } if(!masked_array[f] && vm.coil_residues.includes(f) && topviewer.pluginInstance.domainTypes[index].data[f]) {
            topviewer.pluginInstance.domainTypes[index].data[f].color = "rgb(0,0,0)";
            topviewer.pluginInstance.domainTypes[index].data[f].tooltipMsg = "NaN";
        }                        
        f++;
    }
}
function cleanCustomMap(checked_customMap){
    if (checked_customMap){return;}
    var topviewer = document.getElementById("PdbeTopViewer");
    topviewer.pluginInstance.domainTypes = topviewer.pluginInstance.domainTypes.filter(obj => {return obj.label !== "CustomData"})
    window.coilsOutOfCustom = null;
    //window.custom_prop = null;
    vm.csv_data = null;
}
function handleCustomMappingData(){
    const readFile = function (fileInput) {
        var reader = new FileReader();
        reader.onload = function () {
            vm.csv_data = reader.result;
        };
        reader.readAsBinaryString(fileInput);
    };
    readFile(vm.$refs.custom_csv_file.files[0]);

}
function downloadCSVData() {

    let csv = generateCSVstring(mapped_aa_properties);

    let anchor = document.createElement('a');
    anchor.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv);
    anchor.target = '_blank';
    anchor.download = 'rv3data.csv';
    anchor.click();

}
function handlePropensities(checked_propensities){
    if (checked_propensities){
        console.log("Checked")
    }else{
        console.log("UnChecked")
    }
    
}
function cleanFilter(checked_filter, masking_range){
    if (checked_filter){return;}
    if (masking_range == null){return;}
    window.masked_array = [];
    vm.masking_range = null;
    var topviewer = document.getElementById("PdbeTopViewer");
    topviewer.pluginInstance.getAnnotationFromRibovision(mapped_aa_properties);
    if(window.custom_prop) {
        topviewer.pluginInstance.getAnnotationFromRibovision(window.custom_prop);
    }
    var selectedIndex = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox').selectedIndex;
    topviewer.pluginInstance.updateTheme(topviewer.pluginInstance.domainTypes[selectedIndex].data); 
    window.viewerInstance.visual.select({data: selectSections_RV1.get(topviewer.pluginInstance.domainTypes[selectedIndex].label), nonSelectedColor: {r:255,g:255,b:255}});
}
function cleanSelection(checked_selection, filter_range){
    if (checked_selection){return;}
    if (filter_range == null){return;}
    vm.filter_range = null;
    window.filterRange = "-10000,10000";
    viewerInstance.visual.update({
        customData: {
            url: `https://www.ebi.ac.uk/pdbe/coordinates/${window.pdblower}/chains?entityId=${topviewer.entityId}&encoding=bcif`,
            format: 'cif',
            binary:true },
        assemblyId: '1',
        subscribeEvents: true});
}
function getCookie(name) {
    var cookieValue = null;
    if (document.cookie && document.cookie !== '') {
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = jQuery.trim(cookies[i]);
            // Does this cookie string begin with the name we want?
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
                break;
            }
        }
    }
    return cookieValue;
}
var csrftoken = getCookie('csrftoken');    function csrfSafeMethod(method) {
    // these HTTP methods do not require CSRF protection
    return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
}
$.ajaxSetup({
    beforeSend: function(xhr, settings) {
        if (!csrfSafeMethod(settings.type) && !this.crossDomain) {
            xhr.setRequestHeader("X-CSRFToken", csrftoken);
        }
    }
});

function ajax(url, optional_data='') {
    if (optional_data != ''){
        //var el = document.getElementsByName("csrfmiddlewaretoken");
        //csrf_value = Cookies.get('csrftoken');
        //csrf_value = el[0].getAttribute("value");
        return new Promise((resolve, reject) => {
            $.ajax({
                url: url,
                type: 'POST',
                dataType: "json",
                data: optional_data,
                headers: {'X-CSRFToken': csrftoken},
                success: function(data) {
                    resolve(data)
                },
                error: function(error) {
                    console.log(`Error ${error}`);
                    reject(error)
                }
            })
        })
    }else{
        return new Promise((resolve, reject) => {
            $.ajax({
                url: url,
                type: 'GET',
                dataType: "json",
                success: function(data) {
                    resolve(data)
                },
                error: function(error) {
                    console.log(`Error ${error}`);
                    reject(error)
                }
            })
        })
    }
}

var pushChainData = function(temp_arr, chain_listI){
    temp_arr.push({
        text: chain_listI["molecule_name"][0],
        value: chain_listI["in_chains"][0],
        sequence: chain_listI["sequence"],
        entityID: chain_listI["entity_id"],
        startIndex: chain_listI.source[0].mappings[0].start.residue_number
    })
    return temp_arr;
}

var filterAvailablePolymers = function(chain_list, aln_id, vueObj) {
    let temp_arr = [];
    let url = `/desire-api/alignments/${aln_id}/?format=json`;
    ajax(url).then( aln_data => {
        for (let i = 0; i < chain_list.length; i++) {
            let chain_listI = chain_list[i]
            if (chain_listI["molecule_type"].toLowerCase() == "bound") {continue;}
            if (chain_listI["molecule_type"].toLowerCase() == "water") {continue;}
            for (let ix =0; ix < aln_data["polymers"].length; ix++){
                if (aln_data["polymers"][ix]["genedescription"].trim() == chain_list[i]["molecule_name"][0]){
                    temp_arr = pushChainData(temp_arr, chain_listI);
                }
            }
        }
    // console.log("___" + temp_arr[temp_arr.length - 1]["sequence"] + "___");
    let chain_options = Array.from(new Set(temp_arr.map(JSON.stringify))).map(JSON.parse);
    if (chain_options.length === 0) {
        chain_options.push({text: "Couldn't find polymers from this structure!", value: null})
    }
    vueObj.chains = chain_options;
    });
}

var create_deleted_element = function (parent_id, child_id, child_text) {
    const parent = document.getElementById(parent_id);
    const child_elt = document.createElement("div");
    const childText = document.createTextNode(child_text);
    child_elt.setAttribute("id", child_id);
    child_elt.setAttribute("id", child_id);
    child_elt.appendChild(childText);
    parent.appendChild(child_elt);
}

var cleanupOnNewAlignment = function (vueObj, aln_text='') {
    const menu_item = document.querySelector(".smenubar");
    const aln_item = document.getElementById("alnDiv");
    const topview_item = document.getElementById("topview");
    const molstar_item = document.getElementById("pdbeMolstarView");
    const pdb_input = document.getElementById("pdb_input");
    if (menu_item) {menu_item.remove();}
    if (aln_text != ''){
        vueObj.custom_aln_twc_flag == null;
        window.mapped_aa_properties == null;
        if (pdb_input) {
            if (pdb_input.getAttribute("value") != ""){vueObj.pdbid = null;}
        }
        if (vueObj.chains) {vueObj.chains = null;}
        if (vueObj.aln_meta_data) {vueObj.aln_meta_data = null;}
        if (vueObj.fasta_data) {vueObj.fasta_data = null;}
        if (vueObj.fastaSeqNames) {vueObj.fastaSeqNames = null;}
        if (vueObj.frequency_data) {vueObj.frequency_data = null;}
        if (vueObj.topology_loaded) {vueObj.topology_loaded = 'False';}
        if (aln_item) {aln_item.remove(); create_deleted_element("alnif", "alnDiv", aln_text)}
    }
    window.ajaxRun = false;
    if (window.masked_array.length > 0) {window.masked_array = [];}
    if (vueObj.masking_range) {vueObj.masking_range = null;}
    //if (vueObj.chainid) {vueObj.chainid = null;}
    if (vueObj.checked_filter) {vueObj.checked_filter = false;}
    if (vueObj.checked_customMap) {vueObj.checked_customMap = false;}
    if (vueObj.csv_data) {vueObj.csv_data = null;}
    if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "Select new chain!")}
    if (molstar_item) {molstar_item.remove(); create_deleted_element("molif", "pdbeMolstarView", "Select new structure!")}
}

var loadParaOptions = function (action, callback, vm) {
    if (action === "LOAD_ROOT_OPTIONS"){
        ajax('/alignments/showStrucTaxonomy').then(data =>{
            data.isDisabled = true,
            vm.options = [data];
            callback();
        }).catch(error => {
            console.log(error)
        })
    }
}

var loadParaAlns = function (value, vm) {
    vm.alignments = null;
    ajax('/alignments/fold-api/'+value).then(data=>{
        var fpa = data["Folds to polymers to alignments"]
        var fpa_viz = [];
        Object.keys(fpa).forEach(fkey => {
            Object.keys(fpa[fkey]).forEach(pkey => {
                fpa[fkey][pkey].forEach(function (akey){
                    fpa_viz.push({
                        text:  'Alignment '.concat(akey[1],'; fold ',fkey),
                        value: fkey.concat(',',akey)
                    });
                });
            });
        });
        var temp_arr = fpa_viz
        fpa_viz = Array.from(new Set(temp_arr.map(JSON.stringify))).map(JSON.parse);
        vm.alignments = fpa_viz
    });
}

var calculateFrequencyData = function (frequencies){
    const multiplyvector = function (a,b){
        return a.map((e,i) => e * b[i]);
    }
    let aaPropertiesData = new Map([
                            ["Charge",[0,0,-1,-1,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0]],
                            ["Hydropathy",[1.8,2.5,-3.5,-3.5,2.8,-0.4,-3.2,4.5,-3.9,3.8,1.9,-3.5,-1.6,-3.5,-4.5,-0.8,-0.7,4.2,-0.9,-1.3]],
                            ["Hydrophobicity",[0.02,0.77,-1.04,-1.14,1.35,-0.80,0.26,1.81,-0.41,1.14,1,-0.77,-0.09,-1.10,-0.42,-0.97,-0.77,1.13,1.71,1.11]],
                            ["Polarity",[0,1.48,49.7,49.9,0.35,0,51.6,0.13,49.5,0.13,1.43,3.38,1.58,3.53,52,1.67,1.66,0.13,2.1,1.61]],
                            ["Mutability",[100,44,86,77,51,50,91,103,72,54,93,104,58,84,83,117,107,98,25,50]],
                            ["Shannon entropy",[0.000000000000001,4.321928094887363]],
                            ["TwinCons",[-2.935,12.065]]
                        ]);
    let aaColorData = new Map([
                            ["Charge",[Blues, Reds]],
                            ["Hydropathy",[Blues, Reds]],
                            ["Hydrophobicity",[Reds, Blues]],
                            ["Polarity",[viridis]],
                            ["Mutability",[viridis]],
                            ["Shannon entropy",[plasma]],
                            ["TwinCons",[RdPu, YlGn]],
                        ]);
    window.aaColorData = aaColorData;
    window.aaPropertyConstants = aaPropertiesData;
    window.selectSections_RV1 = new Map();
    let outPropertyPosition = new Map();
    aaPropertiesData.forEach(function (data, property_name){
        if (property_name == "TwinCons"){return;}
        let const_data = data
        outPropertyPosition.set(property_name, [])
        frequencies.forEach(function (col_frequency) {
            if (property_name == "Shannon entropy"){
                const_data = new Array;
                col_frequency.forEach( function (single_freq){
                    if (single_freq == 0){
                        const_data.push(0)
                    }else{
                        const_data.push(Math.log2(single_freq)*-1)
                    }
                });
            }
            outPropertyPosition.get(property_name).push(multiplyvector(const_data, col_frequency));
        });
    });
    return outPropertyPosition;
}

var mapAAProps = function (aa_properties, mapping){
    let outPropertyMappedPosition = new Map();
    aa_properties.forEach(function (data, property_name){
        outPropertyMappedPosition.set(property_name, [])
        data.forEach(function (data, aln_ix) {
            let mappedI0 = mapping[aln_ix+1];
            if (mappedI0) {
                outPropertyMappedPosition.get(property_name).push([mappedI0, Number(math.sum(data).toFixed(2))]);
            }
        });
    });
    return outPropertyMappedPosition;
}

var filterCoilResidues = function (coil_data){
    const range = (start, stop, step) => Array.from({ length: (stop - start) / step + 1}, (_, i) => start + (i * step));
    let coilResidues = [];
    coil_data.forEach(function (coilRange){
        if (coilRange.start < coilRange.stop){
            coilResidues.push(range(coilRange.start, coilRange.stop, 1))
        }
    })
    return coilResidues.flat()
}

var generateCSVstring = function (mapped_data){
    let properties = Array.from(mapped_data.keys());
    let csv = 'Index,'
    csv += properties.join(',');
    csv += '\n';
    let csv_ix = [];
    
    mapped_data.get(properties[0]).forEach((datapoint) =>{
        csv_ix.push([datapoint[0]]);
    })

    properties.forEach((prop) => {
        let ix = 0;
        mapped_data.get(prop).forEach((datapoint) =>{
            csv_ix[ix].push(datapoint[1]);
            ix += 1;
        })
    })

    csv_ix.forEach((row) => {
        csv += row.join(',');
        csv += '\n';
    })

    return csv;
}

var masked_array = [];
Vue.component('treeselect', VueTreeselect.Treeselect, )

var vm = new Vue({
    el: '#phylo_tree_dropdown',
    delimiters: ['[[', ']]'],
    data: {
        tax_id: null,
        alnobj: null,
        options: null,
        alignments: null,
        pdbid: null,
        chains: null,
        chainid: null,
        aln_meta_data: null,
        fasta_data: null,
        fastaSeqNames: null,
        hide_chains: null,
        type_tree: "orth",
        aa_properties: null,
        structure_mapping: null,
        file: null,
        custom_aln_twc_flag: null,
        topology_loaded: 'False',
        twc_loaded: false,
        masking_range: null,
        filter_range: null,
        correct_mask: false,
        coil_residues: null,
        checked_filter: false,
        checked_selection: false,
        checked_customMap: false,
        csv_data: null,
        checked_propensities: false,
        helix_residues: null,
    },
    watch: {
        csv_data: function (csv_data) {
            var topviewer = document.getElementById("PdbeTopViewer");
            var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
            if (csv_data == null){
                selectBoxEle.removeChild(selectBoxEle.childNodes[selectBoxEle.options.length-1]);
                topviewer.pluginInstance.resetDisplay();
                return;
            }
            let custom_data = csv_data.split('\n').map(function(e){
                return e.split(',').map(Number);
            })
            if (custom_data[custom_data.length-1] == 0){custom_data.splice(-1,1)}
            if (topviewer != null && topviewer.pluginInstance.domainTypes != undefined){
                let vals = custom_data.map(function(v){ return v[1] });
                let indexes = custom_data.map(function(v){ return v[0] });
                window.aaColorData.set("CustomData", [viridis]);
                window.aaPropertyConstants.set("CustomData", [Math.min(...vals), Math.max(...vals)]);
                let coilsOutOfCustom = this.coil_residues.filter(value => !indexes.includes(value));
                window.coilsOutOfCustom = coilsOutOfCustom;
                var custom_prop = new Map();
                custom_prop.set("CustomData", custom_data);
                topviewer.pluginInstance.getAnnotationFromRibovision(custom_prop);
                window.custom_prop = custom_prop;
                var custom_option = document.createElement("option");
                custom_option.setAttribute("value", selectBoxEle.options.length);
                custom_option.appendChild(document.createTextNode("Custom Data"));
                selectBoxEle.appendChild(custom_option);
                if(vm.correct_mask == 'True') {
                    var j = topviewer.pluginInstance.domainTypes.length-1;
                    colorResidue(j, window.masked_array);
                }
            }
        },
    },
    methods: {
        handleFileUpload(){
            this.file = this.$refs.custom_aln_file.files[0];
            if (this.tax_id != null){this.tax_id = null;}
        },
        submitCustomAlignment(){
            let formData = new FormData();
            formData.append('custom_aln_file', this.file)
            $.ajax({
                url: '/custom-aln-data',
                data: formData,
                cache: false,
                contentType: false,
                processData: false,
                method: 'POST',
                type: 'POST', // For jQuery < 1.9
                success: function(data){
                    cleanupOnNewAlignment(vm, "Loading alignment...");
                    vm.alnobj = "custom";
                    vm.showAlignment(null, null, "upload");
                },
                error: function(error) {
                    alert(`${error.responseText}`);
                }
            });
        },
        cleanTreeOpts() {
            cleanupOnNewAlignment(vm, "Select new alignment!");
            [this.options, this.tax_id, this.alnobj] = [null, null, null];
        }, loadOptions({ action, callback }) {
            if (this.type_tree == "orth"){
                if (action === "LOAD_CHILDREN_OPTIONS") {
                    action = "";
                    callback();
                //     When they figure out LOAD_CHILDREN_OPTIONS with async search
                //     ajax(`/alignments/showTaxonomy-api/${parentNode.id}`).then(data => {
                //         let fetched_data = [data]
                //         parentNode.children = fetched_data[0].children
                //         callback()
                //     }).catch(error => {
                //         parentNode.children = []
                //         console.log(error)
                //         callback(new Error(`Failed to load options: network error: ${error}`))
                //     })
                };
                if (action === "LOAD_ROOT_OPTIONS") {
                    ajax(`/alignments/showTaxonomy-api/0`).then(data => {
                        data.isDisabled = true;
                        this.options = [data];
                        callback();
                    }).catch(error => {
                        console.log(error)
                    })
                };
                if (action === "LOAD_ROOT_OPTIONS") {
                    ajax('/alignments/showTaxonomy').then(data => {
                        if (this.type_tree == "orth"){
                            this.options = null;
                            data.isDisabled = true;
                            this.options = [data];
                            callback();
                        }
                    }).catch(error => {
                        console.log(error)
                    })
                };
            }
            if (this.type_tree == "para"){
                loadParaOptions(action, callback, vm);
            }
        }, loadData (value, type_tree) {
            if (type_tree == "upload"){this.tax_id = null; return;}
            if (value.length == 0){this.tax_id = null; return;}
            cleanupOnNewAlignment(vm, "Select new alignment!");
            if (this.alnobj != null) {this.alnobj = null;}
            if (type_tree == "orth"){
                this.alignments = null;
                var url = '/desire-api/taxonomic-groups/?format=json&taxgroup_id__in=' + value
                ajax(url).then(data => {
                    if (data["results"].length === 2) {
                        function getObjIntersection(o1, o2) {
                            return Object.keys(o1).filter({}.hasOwnProperty.bind(o2));
                        }
                        var alns_first_tax = Object.fromEntries(data["results"][0]["alignment_ids"]);
                        var alns_second_tax = Object.fromEntries(data["results"][1]["alignment_ids"]);
                        var aln_indexes = getObjIntersection(alns_first_tax, alns_second_tax);
                        var fpa = []
                        aln_indexes.forEach(function(alnk) {
                            fpa.push(Array(Number(alnk), alns_first_tax[alnk]))
                        });
                    } else {
                        var fpa = data["results"][0]["alignment_ids"]
                    }
                    var fpa_viz = [];
                    fpa.forEach(function(fkey) {
                        fpa_viz.push({
                            text: fkey[1],
                            value: fkey[0]
                        });
                    });
                    this.alignments = fpa_viz
                });
            }
            if (type_tree == "para"){
                loadParaAlns (value, vm)
            }
        }, getPDBchains(pdbid, aln_id) {
            if (pdbid.length === 4) {
                if (document.querySelector("pdb-topology-viewer") || document.querySelector("pdbe-molstar")) {cleanupOnNewAlignment(vm);}
                this.chains = null
                this.hide_chains = true
                ajax('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/' + pdbid.toLowerCase())
                    .then(struc_data => {
                        var chain_list = struc_data[pdbid.toLowerCase()];
                        if (this.type_tree == "para") {aln_id = aln_id.split(',')[1]}
                        if (this.type_tree != "upload") {
                            filterAvailablePolymers(chain_list, aln_id, vm);
                            this.hide_chains = null;
                        } else {
                            let chain_options = []
                            for (let i = 0; i < chain_list.length; i++) {
                                let chain_listI = chain_list[i]
                                if (chain_listI["molecule_type"].toLowerCase() == "bound") {continue;}
                                if (chain_listI["molecule_type"].toLowerCase() == "water") {continue;}
                                if (typeof(chain_listI.source[0]) === "undefined") {continue;}
                                chain_options = pushChainData(chain_options, chain_listI);
                            }
                            if (chain_options.length === 0) {
                                chain_options.push({text: "Couldn't find polymers from this structure!", value: null})
                            }
                            vm.chains = chain_options;
                            this.hide_chains = null;
                        }
                    }).catch(error => {
                        alert("Problem with parsing the chains:\n" + error)
                    })
            }
        },
        showAlignment(aln_id, taxid, type_tree) {
            cleanupOnNewAlignment(vm, "Loading alignment...");
            if (type_tree == "orth"){
                var url = `/ortholog-aln-api/${aln_id}/${taxid}`}
            if (type_tree == "para"){
                var url = '/paralog-aln-api/'+aln_id.split(',')[1]}
            if (type_tree == "upload"){
                var url = '/custom-aln-data'}
            ajax(url).then(fasta => {
                if (fasta['TwinCons'] != null){
                    this.custom_aln_twc_flag = fasta['TwinCons']
                }
                vm.fastaSeqNames = fasta['Sequence names'];
                window.aaFreqs = fasta['AA frequencies'];
                var main_elmnt = document.querySelector(".alignment_section");
                window.main_elmnt = main_elmnt;
                let seqsForMSAViewer = parseFastaSeqForMSAViewer(fasta['Alignment']);
                var msaOptions = {
                    sequences: seqsForMSAViewer,
                    colorScheme: "clustal2",
                    height: main_elmnt.offsetHeight * 0.9,
                    width: main_elmnt.offsetWidth * 0.75,
                    tileHeight: 18,
                    tileWidth: 18,
                    overflow: "auto",
                };
                window.msaOptions = msaOptions;
                ReactDOM.render(
                    React.createElement(MyMSA, msaOptions),
                    document.getElementById('alnDiv')
                  );
                this.fasta_data = fasta['Alignment'];
                this.aa_properties = calculateFrequencyData(fasta['AA frequencies']);
            })
        }, showTopologyViewer (pdbid, chainid, fasta){
            window.filterRange = "-10000,10000";
            if (document.querySelector("pdb-topology-viewer") || document.querySelector("pdbe-molstar")) {cleanupOnNewAlignment(vm);}
            if (chainid.length > 1){this.chainid = chainid[0];}
            const topview_item = document.getElementById("topview");
            const molstar_item = document.getElementById("pdbeMolstarView");
            if (topview_item) {topview_item.remove(); create_deleted_element("topif", "topview", "Loading topology viewer and conservation data...")}
            if (molstar_item) {molstar_item.remove(); create_deleted_element("molif", "pdbeMolstarView", "Loading Molstar Component...")}
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
            window.pdblower = pdblower;
            let temp = vm.chains.filter(obj => {
                return obj["value"] == chainid;
            })[0];
            let ebi_sequence = temp["sequence"];
            let startIndex = temp["startIndex"];
            // let ebi_sequence = vm.chains[0]["sequence"];
            ajax('/mapSeqAln/', {fasta, ebi_sequence, startIndex}).then(struct_mapping=>{
                this.structure_mapping = struct_mapping;
                var mapped_aa_properties = mapAAProps(this.aa_properties, struct_mapping);
                if ((this.tax_id != null && this.tax_id.length == 2) || (this.custom_aln_twc_flag != null && this.custom_aln_twc_flag == true) || (this.type_tree == 'para')) {
                    ajax('/twc-api/', {fasta}).then(twcDataUnmapped => {
                        const build_mapped_props = function(mapped_props, twcDataUnmapped, structure_mapping){
                            mapped_props.set("TwinCons", [])
                            for (let i = 0; i < twcDataUnmapped.length; i++) {
                                let mappedI0 = structure_mapping[twcDataUnmapped[i][0]];
                                if (mappedI0) {
                                    mapped_props.get("TwinCons").push([mappedI0, twcDataUnmapped[i][1]]);
                                }
                            }
                            return mapped_props;
                        }
                        var topviewer = document.getElementById("PdbeTopViewer");
                        mapped_aa_properties = build_mapped_props(mapped_aa_properties, twcDataUnmapped, this.structure_mapping);
                        window.mapped_aa_properties = mapped_aa_properties;
                        if (topviewer != null && topviewer.pluginInstance.domainTypes != undefined){
                            var empty_props = new Map();
                            let twc_props = build_mapped_props(empty_props, twcDataUnmapped, this.structure_mapping);
                            topviewer.pluginInstance.getAnnotationFromRibovision(twc_props);
                            var selectBoxEle = topviewer.pluginInstance.targetEle.querySelector('.menuSelectbox');
                            var twc_option = document.createElement("option");
                            twc_option.setAttribute("value", selectBoxEle.options.length);
                            twc_option.appendChild(document.createTextNode("TwinCons"));
                            selectBoxEle.appendChild(twc_option);
                        }
                    })
                }
                window.mapped_aa_properties = mapped_aa_properties;
                var topology_url = `https://www.ebi.ac.uk/pdbe/api/topology/entry/${pdblower}/chain/${chainid}`
                ajax(topology_url).then(data => {
                    var entityid = Object.keys(data[pdblower])[0];
                    vm.coil_residues = filterCoilResidues(data[pdblower][entityid][chainid]["coils"])
                    vm.helix_residues = filterCoilResidues(data[pdblower][entityid][chainid]["helices"])
                    var mapping = [];
                    var range_string = minIndex.concat("-").concat(maxIndex);
                    GetRangeMapping(pdbid, chainid, range_string, mapping);
                    let data_string = JSON.stringify(Array.from(mapped_aa_properties.entries())).replaceAll(",[[", ":").replaceAll("]],",";").replaceAll("],[",",");
                    let formatted_data_string = data_string.replaceAll("[","").replaceAll("]","").replaceAll("\"","");
                    var topology_viewer = `<pdb-topology-viewer id="PdbeTopViewer" entry-id=${pdbid} entity-id=${entityid} chain-id=${chainid}	entropy-id=${formatted_data_string} filter-range=${mapping}></pdb-topology-viewer>`
                    document.getElementById('topview').innerHTML = topology_viewer;
                    window.viewerInstanceTop = document.getElementById("PdbeTopViewer");
                    this.topology_loaded = 'True';
                })
            });
        }, showPDBViewer(pdbid, chainid, entityid){
            if (document.querySelector("pdbe-molstar")) {return;}
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
            window.pdblower = pdblower;
            var viewerInstance = new PDBeMolstarPlugin();
            var options = {
                customData: { url: `https://www.ebi.ac.uk/pdbe/coordinates/${pdblower}/chains?entityId=${entityid}&encoding=bcif`, 
                                format: 'cif', 
                                binary:true },
                hideCanvasControls: ["expand", "selection", " animation"],
                assemblyId: '1',
                hideControls: true,
                subscribeEvents: true,
                bgColor: {r:255,g:255,b:255},
            }
            var viewerContainer = document.getElementById('pdbeMolstarView');
            viewerInstance.render(viewerContainer, options);
            window.viewerInstance = viewerInstance;

            document.addEventListener('PDB.topologyViewer.click', (e) => {
                var molstar= viewerInstance;                            
                var chainId=e.eventData.chainId;
                var entityId=e.eventData.entityId;
                var residueNumber=e.eventData.residueNumber;
                var types=e.eventData.type;                            
                molstar.visual.select({
                    data:[
                        {
                            entity_id:entityId,
                            start_residue_number:residueNumber,
                            end_residue_number:residueNumber,
                            color:{r:20, y:100, b:200},
                            focus:false
                        },
                    ],
                })
            })
            document.addEventListener('PDB.topologyViewer.mouseover', (e) => {
                var molstar= viewerInstance;                            
                var chainId=e.eventData.chainId;
                var entityId=e.eventData.entityId;
                var residueNumber=e.eventData.residueNumber;
                var types=e.eventData.type;
                
                molstar.visual.highlight({
                    data:[
                        {
                            entity_id:entityId,
                            start_residue_number:residueNumber,
                            end_residue_number:residueNumber,
                        },
                    ],
                })
            })
            document.addEventListener('PDB.molstar.mouseover', (e) => {
                var eventData = e.eventData;
                let resi_id = eventData.auth_seq_id;
                if(masked_array && masked_array[resi_id] == false) {
                    viewerInstance.plugin.behaviors.interaction.hover._value.current.loci.kind = "empty-loci"
                }
            });
        }         
    }
});

  (function() {
    var mousePos;
    document.onmousemove = handleMouseMove;
    function handleMouseMove(event) {
        var eventDoc, doc, body;
        event = event || window.event; // IE-ism
        // If pageX/Y aren't available and clientX/Y are,
        // calculate pageX/Y - logic taken from jQuery.
        // (This is to support old IE)
        if (event.pageX == null && event.clientX != null) {
            eventDoc = (event.target && event.target.ownerDocument) || document;
            doc = eventDoc.documentElement;
            body = eventDoc.body;
            event.pageX = event.clientX +
              (doc && doc.scrollLeft || body && body.scrollLeft || 0) -
              (doc && doc.clientLeft || body && body.clientLeft || 0);
            event.pageY = event.clientY +
              (doc && doc.scrollTop  || body && body.scrollTop  || 0) -
              (doc && doc.clientTop  || body && body.clientTop  || 0 );
        }
        mousePos = {
          x: event.pageX,
          y: event.pageY
      };
      window.mousePos = mousePos;
    }
  })();

var Tooltip = function (props) {
    const { style, children, ...otherProps } = props;
    const containerStyle = {
      display: "inline-block"
    };
    const tooltipStyle = {
      position: "relative",
      width: "160px"
    };
    const textStyle = {
      color: "#fff",
      fontSize: "14px",
      lineHeight: 1.2,
      textAlign: "center",
      backgroundColor: "#000",
      borderRadius: "3px",
      padding: "7px"
    };
    return (
      <div id="tooltip" style={{ ...containerStyle, ...style }} {...otherProps}>
        <div style={tooltipStyle}>
          <div style={textStyle}>{children}</div>
        </div>
      </div>
    );
  }
  Tooltip.defaultProps = {
    style: {},
  };

window.ajaxRun = false;
function MyMSA() {
    class SimpleTooltip extends React.Component {
        state = { 
            tileWidth: 18,
            tileHeight: 18,
            aaPos: 0,
            seqPos: 0,
            width: main_elmnt.offsetWidth * 0.7,
            height: main_elmnt.offsetHeight * 0.9,
            highlight: null 
        };
        handleResize = () => {
            this.setState({
                width: main_elmnt.offsetWidth * 0.7,
                height: main_elmnt.offsetHeight * 0.9
            });
            //var style = document.querySelector('[data="rv3_style"]');
            //style.innerHTML = ".slider::-webkit-slider-thumb { width: "+main_elmnt.offsetWidth*0.05+"px}"
        };
        componentDidMount() {
            window.addEventListener("resize", this.handleResize);
            //var style = document.querySelector('[data="rv3_style"]');
            //style.innerHTML = ".slider::-webkit-slider-thumb { width: "+main_elmnt.offsetWidth*0.05+"px}"
        };
        componentWillUnmount() {
            window.removeEventListener("resize", this.handleResize);
        };
        onResidueMouseEnter = e => {
            if (vm.topology_loaded == 'True'){
                let resiPos = vm.structure_mapping[e.position];
                if (resiPos !== undefined){
                    viewerInstanceTop.pluginInstance.highlight(resiPos, resiPos);
                    viewerInstance.visual.highlight({
                        data:[{
                                entity_id:vm.entityId,
                                start_residue_number:resiPos,
                                end_residue_number:resiPos,
                            },],
                    });
                }
            }
            if (!window.ajaxRun){
                window.ajaxRun = true;
                if (e.position !== undefined){
                    registerHoverResiData(e, this);
                }
            } else {
                this.setState({ fold: undefined, phase: undefined });
            }
         };
        onResidueMouseLeave = e => {
            if (vm.topology_loaded == 'True'){
                viewerInstanceTop.pluginInstance.clearHighlight();
                viewerInstance.visual.clearHighlight();
            }
            this.setState({ fold: undefined, phase: undefined });
        };
        highlightRegion = () => {
            const highlight = {
                sequences: {
                  from: 0,
                  to: 2
                },
                residues: {
                  from: 2,
                  to: 13
                }
              };
            this.setState({ highlight });
        };
        removeHighlightRegion = () => {
            this.setState({ highlight: null });
        };                
        render() {
            const xPos = this.state.tileWidth * (this.state.aaPos - 1);
            const yPos = this.state.tileHeight * (this.state.seqPos - 1);
            const maxXpos = window.aaFreqs.length - Math.round(((main_elmnt.offsetWidth * 0.7)/this.state.tileWidth))+2;
            const maxYpos = vm.fastaSeqNames.length - Math.round(((main_elmnt.offsetHeight * 0.9)/this.state.tileHeight))+2;
            return (
            <div style={{ display: "flex" }}>
                <div>
                  <input
                    style = {{ 
                        width: main_elmnt.offsetWidth * 0.7+"px",
                        position: "relative",
                        left: main_elmnt.offsetWidth * 0.2+"px"
                    }}
                    type="range"
                    min="0"
                    max={maxXpos}
                    value={this.state.aaPos}
                    onChange={(evt) => this.setState({ aaPos: evt.target.value })}
                    class="slider"
                    id="xPosSlider"
                    />
                <ReactMSAViewer.MSAViewer 
                {...msaOptions}
                ref={(ref) => (this.el = ref)}
                highlight={this.state.highlight}
                width={this.state.width}
                height={this.state.height}
                tileWidth={this.state.tileWidth}
                tileHeight={this.state.tileHeight}
                position={{ xPos, yPos }}
                >
                <div style={{ position: "relative", display: "flex"}}>
                <ReactMSAViewer.Labels style={{
                    width: main_elmnt.offsetWidth * 0.2
                    }}/>
                <div>
                    <ReactMSAViewer.SequenceViewer
                      onResidueMouseEnter={this.onResidueMouseEnter}
                      onResidueMouseLeave={this.onResidueMouseLeave}
                    />
                    {this.state.fold && (
                      <div
                        style={{
                          position: "absolute",
                          opacity: 0.8,
                          ...this.state.tooltipPosition,
                        }}
                      >
                        <Tooltip>
                          Fold: {this.state.fold} <br></br>
                          Phase: {this.state.phase}
                        </Tooltip>
                      </div>
                    )}
                    </div>
                </div>
                {/* <button onClick={() => this.highlightRegion()}>
                Highlight Region {" "}
              </button> */}
                </ReactMSAViewer.MSAViewer>
            </div>
            <input
                style={{ 
                    width: main_elmnt.offsetHeight*0.9+"px",
                }}
                type="range"
                min="0"
                max={maxYpos}
                value={this.state.seqPos}
                onChange={(evt) => this.setState({ seqPos: evt.target.value })}
                class="slider"
                id="yPosSlider"
                />
            </div>
            );
        }
    }
return <SimpleTooltip />;
};

var absolutePosition = function (el) {
    var
        found,
        left = 0,
        top = 0,
        width = 0,
        height = 0,
        offsetBase = absolutePosition.offsetBase;
    if (!offsetBase && document.body) {
        offsetBase = absolutePosition.offsetBase = document.createElement('div');
        offsetBase.style.cssText = 'position:absolute;left:0;top:0';
        document.body.appendChild(offsetBase);
    }
    if (el && el.ownerDocument === document && 'getBoundingClientRect' in el && offsetBase) {
        var boundingRect = el.getBoundingClientRect();
        var baseRect = offsetBase.getBoundingClientRect();
        found = true;
        left = boundingRect.left - baseRect.left;
        top = boundingRect.top - baseRect.top;
        width = boundingRect.right - boundingRect.left;
        height = boundingRect.bottom - boundingRect.top;
    }
    return {
        found: found,
        left: left,
        top: top,
        width: width,
        height: height,
        right: left + width,
        bottom: top + height
    };
}

var parseFastaSeqForMSAViewer = function (fasta){
    let outSeqs = [];
    let arrayFasta = fasta.split('\n').slice(0, -1);
    arrayFasta.map(function(element, index) {
        if (index % 2 == 0){
            let seqName = element.replaceAll('_', ' ').replaceAll('>', '');
            let seqObj = {'name': seqName, 'sequence': arrayFasta[index+1]}
            outSeqs.push(seqObj);
        }
    });
    return outSeqs;
}

var registerHoverResiData = function (e, tooltipObj){
    const strainQuery = '&res__poldata__strain__strain=';
    var url = `/desire-api/residue-alignment/?format=json&aln_pos=${String(Number(e.position) + 1)}&aln=${vm.alnobj.id}${strainQuery}${vm.fastaSeqNames[Number(e.i)]}`
    ajax(url).then(alnpos_data => {
        var alnViewCanvasEle = document.querySelector("#alnDiv canvas:nth-of-type(1)");
        var alnViewLabelsEle = document.querySelector("#alnDiv div:nth-of-type(2)");
        let boundLabelBox = alnViewLabelsEle.getBoundingClientRect();
        let boundingBox = absolutePosition(alnViewCanvasEle);
        let relativeBox = alnViewCanvasEle.getBoundingClientRect();
        if (alnpos_data.count != 0){
        ajax('/resi-api/' + alnpos_data["results"][0]["res"].split("/")[5]).then(resiData => {
            if (boundingBox.top < mousePos.y && mousePos.y < boundingBox.bottom && boundingBox.left < mousePos.x && mousePos.x < boundingBox.right){
                let tooltipPosition = {
                  top: mousePos.y-boundingBox.top+5 +"px",
                  left: mousePos.x-relativeBox.left+boundLabelBox.right-boundLabelBox.left+5 +"px",
                };
                if (resiData["Structural fold"][0] !== undefined && resiData["Associated data"][0] !== undefined){
                    tooltipObj.setState({
                    fold: resiData["Structural fold"][0][1],
                    phase: resiData["Associated data"][0][1],
                    tooltipPosition,
                  });
                }else{
                    tooltipObj.setState({
                    fold: 'NA',
                    phase: 'NA',
                    tooltipPosition,
                  });
                }
           }
           window.ajaxRun = false;
        });
        }else{
            if (boundingBox.top < mousePos.y && mousePos.y < boundingBox.bottom && boundingBox.left < mousePos.x && mousePos.x < boundingBox.right){
                let tooltipPosition = {
                    top: mousePos.y-boundingBox.top+5 +"px",
                    left: mousePos.x-relativeBox.left+boundLabelBox.right-boundLabelBox.left+5 +"px",
                };
                window.ajaxRun = false;
                tooltipObj.setState({
                    fold: 'NA',
                    phase: 'NA',
                    tooltipPosition,
                });
            }
        }
    }).catch(error => {
       console.log(error);
    })
    return true;
}