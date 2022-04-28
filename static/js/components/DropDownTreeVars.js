export function initialState (){
    return {
        tax_id: null,
        alnobj: null,
        options: null,
        alignments: null,
        pdbs: [
            {id: "4v9d", name: "4V9D E. coli"},
            {id: "4v6u", name: "4V6U P. furiosus"},
            {id: "4v6x", name: "4V6X H. sapiens"},
            ],
        availColorschemes: [
            "buried","cinema","clustal","clustal2","helix","lesk","mae","strand","taylor","turn","zappo",
            ],
        available_properties: [
            {Name:"Charge", url:"static/alignments/svg/Charge.svg"},
            {Name:"Hydropathy", url:"static/alignments/svg/Hydropathy.svg"},
            {Name:"Hydrophobicity", url:"static/alignments/svg/Hydrophobicity.svg"},
            {Name:"Polarity", url:"static/alignments/svg/Polarity1.svg"},
            {Name:"Mutability", url:"static/alignments/svg/Mutability.svg"},
            {Name:"Shannon entropy", url:"static/alignments/svg/Shannon.svg"},
            ],
        domain_list: null,
        selected_domain: [],
        pdbid: null,
        chains: null,
        chainid: [],
        unfilteredChains: null,
        entityID: null,
        fasta_data: null,
        fastaSeqNames: null,
        colorScheme: 'clustal2',
        colorSchemeData: null,
        msavWillMount: null,
        aaPos: 0, seqPos: 0,
        hide_chains: null,
        type_tree: "orth",
        aa_properties: null,
        structure_mapping: null,
        poor_structure_map: null,
        file: null,
        custom_aln_twc_flag: null,
        unmappedTWCdata: null,
        topology_loaded: false,
        twc_loaded: false,
        masking_range: null,
        filter_range: null,
        correct_mask: null,
        domain_or_selection: null,
        checked_domain: false,
        checked_filter: false,
        checked_selection: false,
        checked_customMap: false,
        selected_property: null,
        selected_property_clone: null,
        protein_types: [],
        protein_type_obj: null,
        csv_data: null,
        custom_headers: [],
        raiseCustomCSVWarn: null,
        checked_propensities: false,
        checked_permutation: false,
        all_residues: null,
        coil_residues: null,
        helix_residues: null,
        strand_residues: null,
        substructures: null,
        property: null,
        uploadSession: false,
        fetchingPDBwithCustomAln: null,
        blastPDBresult: [],
        blastMAPresult: null,
        guideOff: true,
        postedPDBEntities: false,
        pdbStart: null,
        pdbEnd: null,
        pdbSeq: null,
        downloadAlignmentOpt: null,
        downloadMapDataOpt: null,
        freqCSV: null,
        customPDBsuccess: null,
        customPDBid: null,
        PDBparsing: false,
        checkedRNA: false,
        schemesMgr: null,
        didCDHit_truncate: null,
        cdhitSelectedOpt: null,
        cdHITReport: null,
        cdhitOpts: [
            {Name:'Download cdhit report', value:'download'},
            {Name:'Reload original alignment', value:'untrunc'}
        ],
    }
}