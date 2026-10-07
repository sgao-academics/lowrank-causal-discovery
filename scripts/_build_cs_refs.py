# -*- coding: utf-8 -*-
"""Emit cs_refs.tex -- the reference list for cs_main.tex.

Current Science numbers references by order of first citation in the text, so
the order is *derived from the manuscript* rather than typed: this script scans
cs_main.tex for \\cite{...}, keeps first-appearance order, and writes a
thebibliography in exactly that order.

House style: Surname, A. B., Surname2, C. D. and Surname3, E. F., Title in
sentence case. Journal abbreviation, year, volume, pages. doi:10.xxxx/yyyy

Every DOI below was resolved and matched against the publisher / DataCite
metadata on 2026-10-07; non-journal items (conference papers, preprints) carry
the arXiv DOI, which is the registered DOI for that work.

Run:  python _build_cs_refs.py
"""
import os
import re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX = os.path.join(HERE, "cs_main.tex")
OUT = os.path.join(HERE, "cs_refs.tex")

E = {
"bray2024globocan":
 "Bray, F., Laversanne, M., Sung, H., Ferlay, J., Siegel, R. L., Soerjomataram, I. and "
 "Jemal, A., Global cancer statistics 2022: GLOBOCAN estimates of incidence and mortality "
 "worldwide for 36 cancers in 185 countries. \\textit{CA Cancer J. Clin.}, 2024, "
 "\\textbf{74}, 229--263. doi:10.3322/caac.21834",
"misra2008gingivobuccal":
 "Misra, S., Chaturvedi, A. and Misra, N. C., Management of gingivobuccal complex cancer. "
 "\\textit{Ann. R. Coll. Surg. Engl.}, 2008, \\textbf{90}, 546--553. "
 "doi:10.1308/003588408X301136",
"singh2025areca":
 "Singh, A. G. and Chaturvedi, P., Areca nut and oral cancer. \\textit{Oral Dis.}, 2025, "
 "\\textbf{31}, 1467--1472. doi:10.1111/odi.14943",
"ambatipudi2012driver":
 "Ambatipudi, S., Gerstung, M., Pandey, M., Samant, T., Patil, A., Kane, S., Desai, R. S., "
 "Sch\\\"affer, A. A., Beerenwinkel, N. and Mahimkar, M. B., Genome-wide expression and copy "
 "number analysis identifies driver genes in gingivobuccal cancers. \\textit{Genes "
 "Chromosomes Cancer}, 2012, \\textbf{51}, 161--173. doi:10.1002/gcc.20940",
"pansare2019gingivobuccal":
 "Pansare, K., Gardi, N., Kamat, S., Dange, P., Previn, R., Maheshwari, N., Gupta, S., "
 "Chavan, S., Kane, S. and Mahimkar, M., Establishment and genomic characterization of "
 "gingivobuccal carcinoma cell lines with smokeless tobacco associated genetic alterations "
 "and oncogenic PIK3CA mutation. \\textit{Sci. Rep.}, 2019, \\textbf{9}, 8272. "
 "doi:10.1038/s41598-019-44143-0",
"bhosale2017progression":
 "Bhosale, P. G., Cristea, S., Ambatipudi, S., Desai, R. S., Kumar, R., Patil, A., Kane, S., "
 "Borges, A. M., Sch\\\"affer, A. A., Beerenwinkel, N. and Mahimkar, M. B., Chromosomal "
 "alterations and gene expression changes associated with the progression of leukoplakia to "
 "advanced gingivobuccal cancer. \\textit{Transl. Oncol.}, 2017, \\textbf{10}, 396--409. "
 "doi:10.1016/j.tranon.2017.03.008",
"inchanalkar2023methylation":
 "Inchanalkar, M., Srivatsa, S., Ambatipudi, S., Bhosale, P. G., Patil, A., Sch\\\"affer, A. A., "
 "Beerenwinkel, N. and Mahimkar, M. B., Genome-wide DNA methylation profiling of HPV-negative "
 "leukoplakia and gingivobuccal complex cancers. \\textit{Clin. Epigenetics}, 2023, "
 "\\textbf{15}, 93. doi:10.1186/s13148-023-01510-z",
"squarize2006crosstalk":
 "Squarize, C. H., Castilho, R. M., Sriuranpong, V., Pinto, D. S. and Gutkind, J. S., "
 "Molecular cross-talk between the NFkappaB and STAT3 signaling pathways in head and neck "
 "squamous cell carcinoma. \\textit{Neoplasia}, 2006, \\textbf{8}, 733--746. "
 "doi:10.1593/neo.06274",
"monisha2017nfkb":
 "Monisha, J., Roy, N. K., Bordoloi, D., Kumar, A., Golla, R., Kotoky, J., Padmavathi, G. and "
 "Kunnumakkara, A. B., Nuclear factor kappa B: a potential target to persecute head and neck "
 "cancer. \\textit{Curr. Drug Targets}, 2017, \\textbf{18}, 232--253. "
 "doi:10.2174/1389450117666160201112330",
"rao2010proinflammatory":
 "Rao, S. K., Pavicevic, Z., Du, Z., Kim, J. G., Fan, M., Jiao, Y., Rosebush, M., Samant, S., "
 "Gu, W., Pfeffer, L. M. and Nosrat, C. A., Pro-inflammatory genes as biomarkers and "
 "therapeutic targets in oral squamous cell carcinoma. \\textit{J. Biol. Chem.}, 2010, "
 "\\textbf{285}, 32512--32521. doi:10.1074/jbc.M110.150490",
"nawaz2026training":
 "Nawaz, H., Purohit, B. M., Priya, H. and Unnikrishnan, G., Effectiveness of oral cancer "
 "training programs for frontline health workers: a systematic review and meta-analysis. "
 "\\textit{J. Cancer Policy}, 2026, \\textbf{48}, 100741. doi:10.1016/j.jcpo.2026.100741",
"hanahan2022hallmarks":
 "Hanahan, D., Hallmarks of cancer: new dimensions. \\textit{Cancer Discov.}, "
 "2022, \\textbf{12}, 31--46. doi:10.1158/2159-8290.CD-21-1059",
"weinstein2013tcga":
 "Weinstein, J. N. \\textit{et al.}, The Cancer Genome Atlas pan-cancer analysis "
 "project. \\textit{Nat. Genet.}, 2013, \\textbf{45}, 1113--1120. doi:10.1038/ng.2764",
"hoadley2018celloforigin":
 "Hoadley, K. A. \\textit{et al.}, Cell-of-origin patterns dominate the molecular "
 "classification of 10,000 tumors from 33 types of cancer. \\textit{Cell}, 2018, "
 "\\textbf{173}, 291--304. doi:10.1016/j.cell.2018.03.022",
"barabasi2011network":
 "Barab\\'asi, A.-L., Gulbahce, N. and Loscalzo, J., Network medicine: a "
 "network-based approach to human disease. \\textit{Nat. Rev. Genet.}, 2011, "
 "\\textbf{12}, 56--68. doi:10.1038/nrg2918",
"huynh2010inferring":
 "Huynh-Thu, V. A., Irrthum, A., Wehenkel, L. and Geurts, P., Inferring regulatory "
 "networks from expression data using tree-based methods. \\textit{PLoS ONE}, 2010, "
 "\\textbf{5}, e12776. doi:10.1371/journal.pone.0012776",
"moerman2019grnboost2":
 "Moerman, T. \\textit{et al.}, GRNBoost2 and Arboreto: efficient and scalable "
 "inference of gene regulatory networks. \\textit{Bioinformatics}, 2019, "
 "\\textbf{35}, 2159--2161. doi:10.1093/bioinformatics/bty916",
"marbach2012wisdom":
 "Marbach, D. \\textit{et al.}, Wisdom of crowds for robust gene network inference. "
 "\\textit{Nat. Methods}, 2012, \\textbf{9}, 796--804. doi:10.1038/nmeth.2016",
"cs_glyco_oral":
 "Goni, M., Desai, A. K., Kumar, N., Hegde, V. \\textit{et al.}, Significance of "
 "altered glycosyltransferase expression levels in oral cancer. \\textit{Curr. Sci.}, "
 "2022, \\textbf{123}, 52--58. doi:10.18520/cs/v123/i1/52-58",
"cs_izudheen2020":
 "Izudheen, S., Sajan, E. S., George, I., John, J. and Attipetty, C. S., Effect of "
 "community structures in protein--protein interaction network in cancer protein "
 "identification. \\textit{Curr. Sci.}, 2020, \\textbf{118}, 62--69. "
 "doi:10.18520/cs/v118/i1/62-69",
"pearl2009causality":
 "Pearl, J., Causal inference in statistics: an overview. \\textit{Stat. Surv.}, "
 "2009, \\textbf{3}, 96--146. doi:10.1214/09-SS057",
"spirtes2000causation":
 "Spirtes, P., Glymour, C. N. and Scheines, R., \\textit{Causation, Prediction, and "
 "Search}, MIT Press, Cambridge, 2000, 2nd edn. doi:10.7551/mitpress/1754.001.0001",
"heinze2018causal":
 "Heinze-Deml, C., Maathuis, M. H. and Meinshausen, N., Causal structure learning. "
 "\\textit{Annu. Rev. Stat. Appl.}, 2018, \\textbf{5}, 371--391. "
 "doi:10.1146/annurev-statistics-031017-100630",
"runge2019detecting":
 "Runge, J., Nowack, P., Kretschmer, M., Flaxman, S. and Sejdinovic, D., Detecting "
 "and quantifying causal associations in large nonlinear time series datasets. "
 "\\textit{Sci. Adv.}, 2019, \\textbf{5}, eaau4996. doi:10.1126/sciadv.aau4996",
"scholkopf2021toward":
 "Sch\\\"olkopf, B. \\textit{et al.}, Toward causal representation learning. "
 "\\textit{Proc. IEEE}, 2021, \\textbf{109}, 612--634. doi:10.1109/JPROC.2021.3058954",
"shimizu2014lingam":
 "Shimizu, S., LiNGAM: non-Gaussian methods for estimating causal structures. "
 "\\textit{Behaviormetrika}, 2014, \\textbf{41}, 65--98. doi:10.2333/bhmk.41.65",
"friedman2000using":
 "Friedman, N., Linial, M., Nachman, I. and Pe\\'er, D., Using Bayesian networks to "
 "analyze expression data. \\textit{J. Comput. Biol.}, 2000, \\textbf{7}, 601--620. "
 "doi:10.1089/106652700750050961",
"sachs2005causal":
 "Sachs, K., Perez, O., Pe\\'er, D., Lauffenburger, D. A. and Nolan, G. P., Causal "
 "protein-signaling networks derived from multiparameter single-cell data. "
 "\\textit{Science}, 2005, \\textbf{308}, 523--529. doi:10.1126/science.1105809",
"zheng2018dags":
 "Zheng, X., Aragam, B., Ravikumar, P. and Xing, E. P., DAGs with NO TEARS: "
 "continuous optimization for structure learning. In \\textit{Advances in Neural "
 "Information Processing Systems 31}, Curran Associates, Red Hook, 2018. "
 "doi:10.48550/arXiv.1803.01422",
"ng2020golem":
 "Ng, I., Ghassami, A. and Zhang, K., On the role of sparsity and DAG constraints "
 "for learning linear DAGs. In \\textit{Advances in Neural Information Processing "
 "Systems 33}, Curran Associates, Red Hook, 2020. doi:10.48550/arXiv.2006.10201",
"bello2022dagma":
 "Bello, K., Aragam, B. and Ravikumar, P., DAGMA: learning DAGs via M-matrices and "
 "a log-determinant acyclicity characterization. In \\textit{Advances in Neural "
 "Information Processing Systems 35}, Curran Associates, Red Hook, 2022. "
 "doi:10.48550/arXiv.2209.08037",
"lachapelle2019gradient":
 "Lachapelle, S., Brouillard, P., Deleu, T. and Lacoste-Julien, S., Gradient-based "
 "neural DAG learning. In \\textit{International Conference on Learning "
 "Representations}, 2020. doi:10.48550/arXiv.1906.02226",
"yu2019daggnn":
 "Yu, Y., Chen, J., Gao, T. and Yu, M., DAG-GNN: DAG structure learning with graph "
 "neural networks. In \\textit{Proceedings of the 36th International Conference on "
 "Machine Learning}, 2019. doi:10.48550/arXiv.1904.10098",
"vaquerizas2009census":
 "Vaquerizas, J. M., Kummerfeld, S. K., Teichmann, S. A. and Luscombe, N. M., A "
 "census of human transcription factors: function, expression and evolution. "
 "\\textit{Nat. Rev. Genet.}, 2009, \\textbf{10}, 252--263. doi:10.1038/nrg2538",
"lambert2018human":
 "Lambert, S. A. \\textit{et al.}, The human transcription factors. \\textit{Cell}, "
 "2018, \\textbf{172}, 650--665. doi:10.1016/j.cell.2018.01.029",
"kanehisa2021kegg":
 "Kanehisa, M., Furumichi, M., Sato, Y., Ishiguro-Watanabe, M. and Tanabe, M., KEGG: "
 "integrating viruses and cellular organisms. \\textit{Nucleic Acids Res.}, 2021, "
 "\\textbf{49}, D545--D551. doi:10.1093/nar/gkaa970",
"goldman2020ucsc":
 "Goldman, M. J. \\textit{et al.}, Visualizing and interpreting cancer genomics data "
 "via the Xena platform. \\textit{Nat. Biotechnol.}, 2020, \\textbf{38}, 675--678. "
 "doi:10.1038/s41587-020-0546-8",
"tsherniak2017defining":
 "Tsherniak, A. \\textit{et al.}, Defining a cancer dependency map. \\textit{Cell}, "
 "2017, \\textbf{170}, 564--576. doi:10.1016/j.cell.2017.06.010",
"dempster2021chronos":
 "Dempster, J. M. \\textit{et al.}, Chronos: a cell population dynamics model of "
 "CRISPR experiments that improves inference of gene fitness effects. "
 "\\textit{Genome Biol.}, 2021, \\textbf{22}, 343. doi:10.1186/s13059-021-02540-7",
"corsello2020discovering":
 "Corsello, S. M. \\textit{et al.}, Discovering the anticancer potential of "
 "non-oncology drugs by systematic viability profiling. \\textit{Nat. Cancer}, 2020, "
 "\\textbf{1}, 235--248. doi:10.1038/s43018-019-0018-6",
"han2018trrust":
 "Han, H. \\textit{et al.}, TRRUST v2: an expanded reference database of human and "
 "mouse transcriptional regulatory interactions. \\textit{Nucleic Acids Res.}, 2018, "
 "\\textbf{46}, D380--D386. doi:10.1093/nar/gkx1013",
"liberzon2015hallmark":
 "Liberzon, A. \\textit{et al.}, The Molecular Signatures Database (MSigDB) hallmark "
 "gene set collection. \\textit{Cell Syst.}, 2015, \\textbf{1}, 417--425. "
 "doi:10.1016/j.cels.2015.12.004",
"wolf2018scanpy":
 "Wolf, F. A., Angerer, P. and Theis, F. J., SCANPY: large-scale single-cell gene "
 "expression data analysis. \\textit{Genome Biol.}, 2018, \\textbf{19}, 15. "
 "doi:10.1186/s13059-017-1382-0",
"cs_egf_a549":
 "Castro-Aceituno, V., Siddiqi, M. H., Ahn, S., Sathishkumar, N. \\textit{et al.}, "
 "The inhibitory mechanism of compound K on A549 lung cancer cells through EGF "
 "pathway: an \\textit{in silico} and \\textit{in vitro} approach. \\textit{Curr. "
 "Sci.}, 2016, \\textbf{111}, 1071--1077. doi:10.18520/cs/v111/i6/1071-1077",
"wang2018s100a8":
 "Wang, S. \\textit{et al.}, S100A8/A9 in inflammation. \\textit{Front. Immunol.}, "
 "2018, \\textbf{9}, 1298. doi:10.3389/fimmu.2018.01298",
"gottesman2002multidrug":
 "Gottesman, M. M., Fojo, T. and Bates, S. E., Multidrug resistance in cancer: role "
 "of ATP-dependent transporters. \\textit{Nat. Rev. Cancer}, 2002, \\textbf{2}, "
 "48--58. doi:10.1038/nrc706",
}

src = open(TEX, encoding="utf-8").read()
src = re.sub(r"(?m)%.*$", "", src)                 # strip comments
order, seen = [], set()
for m in re.finditer(r"\\cite\{([^}]*)\}", src):
    for k in m.group(1).split(","):
        k = k.strip()
        if k and k not in seen:
            seen.add(k)
            order.append(k)

missing = [k for k in order if k not in E]
assert not missing, "no formatted entry for: %s" % missing
uncited = [k for k in E if k not in seen]
if uncited:
    print("note: declared but never cited ->", uncited)

lines = ["% generated by scripts/_build_cs_refs.py -- do not edit by hand",
         "\\begin{thebibliography}{99}", "\\setlength{\\itemsep}{0pt}",
         "\\small", "\\begin{spacing}{1.0}"]
for i, k in enumerate(order, 1):
    lines.append("\\bibitem{%s} %s" % (k, E[k]))
lines += ["\\end{spacing}", "\\end{thebibliography}", ""]

open(OUT, "w", encoding="utf-8").write("\n".join(lines))
print("wrote %s  (%d references, order = first citation)" % (OUT, len(order)))
for i, k in enumerate(order, 1):
    print("  %2d  %s" % (i, k))
