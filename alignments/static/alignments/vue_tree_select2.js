function ajax(url) {
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
    },
    methods: {
        limiter(e) {
            if (e.length > 2) {
                alert('You can only select two groups!')
                e.pop()
            }
        },
        loadOptions({ action, callback }) {
            if (action === "LOAD_ROOT_OPTIONS") {
                ajax('/alignments/showTaxonomy').then(data => {
                    data.isDisabled = true,
                        this.options = [data];
                    callback();
                }).catch(error => {
                    console.log(error)
                })
            }
        },
        loadData: function(value) {
            this.alignments = null;
            var url = '/desire-api/taxonomic-groups/?format=json&taxgroup_id__in=' + value
            ajax(url)
                .then(data => {
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
        },
        getPDBchains(pdbid) {
            if (pdbid.length === 4) {
                ajax('https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/' + pdbid.toLowerCase())
                    .then(struc_data => {
                        var chain_options = [];
                        for (var i = 0; i < struc_data[pdbid.toLowerCase()].length; i++) {
                            if (struc_data[pdbid.toLowerCase()][i]["molecule_type"] == "Bound") {
                                continue;
                            }
                            if (struc_data[pdbid.toLowerCase()][i]["molecule_type"] == "Water") {
                                continue;
                            }
                            chain_options.push({
                                text: struc_data[pdbid.toLowerCase()][i]["molecule_name"][0],
                                value: struc_data[pdbid.toLowerCase()][i]["in_chains"][0]
                            })
                        }
                        var temp_arr = chain_options
                        chain_options = Array.from(new Set(temp_arr.map(JSON.stringify))).map(JSON.parse);
                        this.chains = chain_options
                    }).catch(error => {
                        alert("No such pdb id: " + pdbid + ".", error)
                    })
            }
        },
        showAlignment(aln_id, taxid) {
            this.aln_meta_data = null;
            var url = '/ortholog-aln-api/' + aln_id + '/' + taxid
            var main_elmnt = document.getElementById("main_elt")
            ajax(url).then(fasta => {
                var opts = {
                    el: document.getElementById("alnDiv"),
                    seqs: msa.io.fasta.parse(fasta[0]),
                    colorscheme: {
                        scheme: "clustal2",
                    },
                    //columns: {
                    //    hidden: fasta[2] // hidden columns
                    //},
                    zoomer: {
                        // general
                        alignmentWidth: main_elmnt.offsetWidth * 0.7,
                        alignmentHeight: main_elmnt.offsetHeight * 0.4,
                        columnWidth: 15,
                        rowHeight: 15,
                        labelNameLength: 300,
                        autoResize: false, // only for the width
                    },
                    conf: {
                        registerMouseHover: false,
                        registerMouseClicks: true,
                    },
                    // smaller menu for JSBin
                    //menu: "small",
                    //bootstrapMenu: true
                };
                var m = new msa.msa(opts);
                m.render();
                m.g.on("residue:click", function(data) {
                    vm.aln_meta_data = null;
                    const strainQuery = '&res__poldata__strain__strain=';
                    var url = `/desire-api/residue-alignment/?format=json&aln_pos=${String(Number(data["rowPos"]) + 1)}&aln=${aln_id}${strainQuery}${fasta[1][Number(data["seqId"])]}`
                    ajax(url).then(alnpos_data => {
                        ajax('/resi-api/' + alnpos_data["results"][0]["res"].split("/")[5]).then(resiData => {
                            vm.aln_meta_data = resiData;
                        });
                    }).catch(error => {
                        vm.aln_meta_data = null;
                        console.log("No residue with alignment position: " + data["rowPos"] + ". In alignment " + aln_id + ". Of species " + fasta[1][Number(data["seqId"])]);
                        console.log(error);
                    })
                });
                m.g.on("residue:mousein", function(data) {
                    console.log(data)
                });
            })
        }, showTopologyViewer (pdbid, chainid, entropy_address){
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
            var topology_url = `https://www.ebi.ac.uk/pdbe/api/topology/entry/${pdblower}/chain/${chainid}`
            ajax(topology_url).then(data => {
                var entityid = Object.keys(data[pdblower])[0];
                var mapping = []
                var range_string = minIndex.concat("-").concat(maxIndex)
                GetRangeMapping(pdbid, chainid, range_string, mapping)
                console.log(mapping)
                var topology_viewer = `<pdb-topology-viewer entry-id=${pdbid} entity-id=${entityid} chain-id=${chainid}	entropy-id=${entropy_address} filter-range=${mapping}></pdb-topology-viewer>`
                document.getElementById('topview').innerHTML = topology_viewer;
            });
        }, 
        

        
        showPDBViewer(pdbid, chainid){
            var minIndex = String(0)
            var maxIndex = String(100000)
            var pdblower = pdbid.toLocaleLowerCase();
            console.log('PDBV');
            var coordinates_url=`https://www.ebi.ac.uk/pdbe/coordinates/${pdblower}/chains?entityId=27&encoding=bcif`;
            var topology_url = `https://www.ebi.ac.uk/pdbe/api/topology/entry/${pdblower}/chain/${chainid}`
            console.log(coordinates_url);

            ajax(topology_url).then (data => {
                var entityid = Object.keys(data[pdblower])[0];
                var coordinates_url=`https://www.ebi.ac.uk/pdbe/coordinates/${pdblower}/chains?${entityid}&encoding=bcif`;
                console.log(entityid);
                var mapping = []
                var range_string = minIndex.concat("-").concat(maxIndex)
                GetRangeMapping(pdbid, chainid, range_string, mapping)
                console.log(mapping)

                
                var PDBMolstar_viewer = `<pdbe-molstar id="PdbeMolstarComponent" molecule-id="1cbs" hide-controls="true" subscribe-events="true" ></pdbe-molstar>`
                document.getElementById('pdbeMolstarView').innerHTML = PDBMolstar_viewer;
                var PdbeMolstarComponent = document.getElementById('PdbeMolstarComponent');
                var viewerInstance2 = PdbeMolstarComponent.viewerInstance;


                viewerInstance2.visual.update({
                    customData: { url: `https://www.ebi.ac.uk/pdbe/coordinates/${pdblower}/chains?entityId=${entityid}&encoding=bcif`, format: 'cif', binary:true },
                          hideCanvasControls: ["expand", "selection", " animation"],
                          assemblyId: '1',                    
                          hideControls: true,                   
                          subscribeEvents: true
                        });

  

                        document.addEventListener('PDB.topologyViewer.click', (e)=>{
                            var pdbeMolstar=document.getElementById("PdbeMolstarComponent")
                            var molstar= pdbeMolstar.viewerInstance;                            
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
                            focus:false},
                            
                            
                            ],
                            
                            })
                            })
                            
                            
                            document.addEventListener('PDB.topologyViewer.mouseover', (e)=>{
                            var pdbeMolstar=document.getElementById("PdbeMolstarComponent")
                            var molstar= pdbeMolstar.viewerInstance;                            
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
                                     


            });
        }


        
    }
 
})


   
   
   
   
   
    