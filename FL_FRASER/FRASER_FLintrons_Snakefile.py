import pandas as pd
import os

path = "/s/project/fraser/fraser2/GTEx_v8/FRASER2_results/minK20_25_minN10/PCA__pc0.1"

tissues = [tissue for tissue in os.listdir(path) if os.path.isdir(os.path.join(path, tissue))] #list of tissues to analyze

rule all:
    input:
        "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_first_introns.csv", 
        "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_last_introns.csv",
        "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_first_introns_all.csv",
        "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_last_introns_all.csv",
        "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_first_introns_withnonmerged.csv",
        "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_last_introns_withnonmerged.csv"
        
rule analyze_tissue:
    input: 
        fraser_results_path = "/s/project/fraser/fraser2/GTfiEx_v8/FRASER2_results/minK20_25_minN10/PCA__pc0.1/{tissue}/optQ__newFilt/delta0.1/results_junction.tsv",
        fds_path = "/s/project/first_last_exon/Data/fds_files/fds_{tissue}.tsv"
    output:
        first_mergedoutput = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_first_introns.csv",
        last_mergedoutput = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_last_introns.csv",
        first_output = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_first_introns_all.csv",
        last_output = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_last_introns_all.csv",
        first_withnonmerged = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_first_introns_withnonmerged.csv",
        last_withnonmerged = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_last_introns_withnonmerged.csv"
    run:
        import pandas as pd
        fraser_results = pd.read_csv(input.fraser_results_path, sep="\t")
        fraser_results = fraser_results.assign(hgncSymbol=fraser_results.hgncSymbol.str.split(';')).explode('hgncSymbol').reset_index(drop=True)
    
        fds = pd.read_csv(input.fds_path, sep="\t")
        #separate entries with more than one gene attributions
        fds_sepsymbols = fds.assign(hgnc_symbol=fds.hgnc_symbol.str.split(';')).explode('hgnc_symbol').reset_index(drop=True)
        fds_sepsymbols.sort_values("start")
        
        first_introns = fds_sepsymbols.groupby("hgnc_symbol").apply(lambda x: x.iloc[0] if x.iloc[0]['strand'] == '+' else x.iloc[-1])
        first_introns.to_csv(output.first_output, index=False)
        first_introns_matched = pd.merge(fraser_results, first_introns, on = ['seqnames', 'start', 'strand'])
        first_introns_matched.to_csv(output.first_mergedoutput, index=False)
        first_introns_aftermerge = pd.merge(fraser_results, first_introns, on=['seqnames', 'start', 'strand'], how='left')
        first_introns_aftermerge.to_csv(output.first_withnonmerged, index=False)
        
        last_introns = fds_sepsymbols.groupby("hgnc_symbol").apply(lambda x: x.iloc[-1] if x.iloc[-1]['strand'] == '+' else x.iloc[0])
        last_introns.to_csv(output.last_output, index=False)
        last_introns_matched = pd.merge(fraser_results, last_introns, on = ['seqnames', 'end', 'strand'])
        last_introns_matched.to_csv(output.last_mergedoutput, index=False)
        last_introns_aftermerge = pd.merge(fraser_results, last_introns, on=['seqnames', 'end', 'strand'], how='left')
        last_introns_aftermerge.to_csv(output.last_withnonmerged, index=False)
        

rule generate_fds:
    output: "/s/project/first_last_exon/Data/fds_files/fds_{tissue}.tsv"
    shell:
        """
        Rscript generate_fds.R {wildcards.tissue}
        """
        
rule merge_results:
    input: 
        first = expand("/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_first_introns.csv", tissue = tissues),
        last = expand("/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_last_introns.csv", tissue = tissues),
        first_all = expand("/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_first_introns_all.csv", tissue = tissues),
        last_all = expand("/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_last_introns_all.csv", tissue = tissues),
        first_withnonmerged = expand("/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_first_introns_withnonmerged.csv", tissue = tissues),
        last_withnonmerged = expand("/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_last_introns_withnonmerged.csv", tissue = tissues)
        
    output:
        first = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_first_introns.csv",
        last = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_last_introns.csv",
        first_all = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_first_introns_all.csv",
        last_all = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_last_introns_all.csv",
        first_withnonmerged = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_first_introns_withnonmerged.csv",
        last_withnonmerged = "/s/project/first_last_exon/Data/FLexons_inFRASERresults/merged_last_introns_withnonmerged.csv"
        
    run:
        dfs_first = []
        dfs_last = []
        dfs_first_all = []
        dfs_last_all = []
        dfs_first_withnonmerged = []
        dfs_last_withnonmerged = []
        for tissue in tissues:
            first = pd.read_csv(f"/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_first_introns.csv")
            last = pd.read_csv(f"/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_last_introns.csv")
            first_all = pd.read_csv(f"/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_first_introns_all.csv")
            last_all = pd.read_csv(f"/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_last_introns_all.csv")
            first_withnonmerged = pd.read_csv(f"/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_first_introns_withnonmerged.csv")
            last_withnonmerged = pd.read_csv(f"/s/project/first_last_exon/Data/FLexons_inFRASERresults/{tissue}_results_last_introns_withnonmerged.csv")
            
            first["tissue"] = tissue
            last["tissue"] = tissue
            first_all["tissue"] = tissue
            last_all["tissue"] = tissue
            first_withnonmerged["tissue"] = tissue
            last_withnonmerged["tissue"] = tissue
            
            dfs_first.append(first)
            dfs_last.append(last)
            dfs_first_all.append(first_all)
            dfs_last_all.append(last_all)
            dfs_first_withnonmerged.append(first_withnonmerged)
            dfs_last_withnonmerged.append(last_withnonmerged)
            
        merged_first = pd.concat(dfs_first)
        merged_last = pd.concat(dfs_last)
        merged_first_all = pd.concat(dfs_first_all)
        merged_last_all = pd.concat(dfs_last_all)
        merged_first_withnonmerged = pd.concat(dfs_first_withnonmerged)
        merged_last_withnonmerged = pd.concat(dfs_last_withnonmerged)
        
        merged_first.to_csv(f"{output.first}", index=False)
        merged_last.to_csv(f"{output.last}", index=False)
        merged_first_all.to_csv(f"{output.first_all}", index=False)
        merged_last_all.to_csv(f"{output.last_all}", index=False)
        merged_first_withnonmerged.to_csv(f"{output.first_withnonmerged}", index=False)
        merged_last_withnonmerged.to_csv(f"{output.last_withnonmerged}", index=False)

