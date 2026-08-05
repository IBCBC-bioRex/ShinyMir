#### Shared palette for ancestor/pathway group colouring (10 colours, wraps cyclically)
GROUP_PALETTE <- c(
  "#FFB3B3", "#B3D1FF", "#B3FFB8", "#FFF3B3",
  "#D4B3FF", "#B3FFE8", "#FFD9B3", "#F3B3FF",
  "#B3B3FF", "#FFDDB3"
)

####queries
query_md <- list(
  query_db_general = "
    SELECT DISTINCT
      A.MIRNA_PREMATURE AS MIRNA_NAME,
      C.DISEASE AS DISEASE,
      B.PUBMED_ID AS PUBMED_ID
    FROM (
      select MIRNA_ID,MIRNA_PREMATURE from MIRNAS A 
    ) A JOIN MIRNAS_DISEASES_ARTICLES B  ON A.MIRNA_ID = B.MIRNA_ID
    JOIN (
      select DISEASE_ID,DISEASE from DISEASES C 
    ) C ON B.DISEASE_ID = C.DISEASE_ID
  ",
  query_general = "
    SELECT DISTINCT
      A.MIRNA_PREMATURE AS MIRNA_NAME,
      C.DISEASE AS DISEASE,
      B.PUBMED_ID AS PUBMED_ID
    FROM (
      select MIRNA_ID,MIRNA_PREMATURE from MIRNAS A where 1=1 {mirna_clause}
    ) A JOIN MIRNAS_DISEASES_ARTICLES B  ON A.MIRNA_ID = B.MIRNA_ID
    JOIN (
      select DISEASE_ID,DISEASE from DISEASES C WHERE 1=1 {disease_clause}
    ) C ON B.DISEASE_ID = C.DISEASE_ID
    {min_assoc_clause}
    WHERE 1=1 {sample_type_clause} {human_clause}
  ",
  query_mirna = "
      SELECT DISTINCT
        A.MIRNA_PREMATURE AS MIRNA_NAME,
        ANY_VALUE(A.MIRBASE_ACC) AS MIRBASE_ACC,
        COUNT(DISTINCT C.DISEASE) AS DISEASE_COUNT,
        GROUP_CONCAT(DISTINCT  C.DISEASE) AS DISEASE,
        COUNT(DISTINCT B.PUBMED_ID) AS PUBMED_COUNT,
        GROUP_CONCAT(DISTINCT  B.PUBMED_ID) AS PUBMED_ID
      FROM (select MIRNA_ID,MIRNA_PREMATURE,MIRBASE_ACC from MIRNAS A where 1=1 {mirna_clause}) A
      JOIN MIRNAS_DISEASES_ARTICLES B ON A.MIRNA_ID = B.MIRNA_ID
      JOIN (select DISEASE_ID,DISEASE from DISEASES C WHERE 1=1 {disease_clause}) C ON B.DISEASE_ID = C.DISEASE_ID
      {min_assoc_clause}
      WHERE 1=1 {sample_type_clause} {human_clause}
      GROUP BY MIRNA_NAME
      ORDER BY DISEASE_COUNT DESC
    ",
  query_disease = "
    SELECT DISTINCT
      C.DISEASE AS DISEASE,
      COUNT(DISTINCT A.MIRNA_PREMATURE) AS MIRNA_COUNT,
      GROUP_CONCAT(DISTINCT  A.MIRNA_PREMATURE) AS MIRNA_NAME,
      COUNT(DISTINCT B.PUBMED_ID) AS PUBMED_COUNT,
      GROUP_CONCAT(DISTINCT  B.PUBMED_ID) AS PUBMED_ID
    FROM (select MIRNA_ID,MIRNA_PREMATURE,MIRBASE_ACC from MIRNAS A where 1=1 {mirna_clause}) A
    JOIN MIRNAS_DISEASES_ARTICLES B ON A.MIRNA_ID = B.MIRNA_ID
    JOIN (select DISEASE_ID,DISEASE from DISEASES C WHERE 1=1 {disease_clause}) C ON B.DISEASE_ID = C.DISEASE_ID
    {min_assoc_clause}
    WHERE 1=1 {sample_type_clause} {human_clause}
    GROUP BY C.DISEASE
    ORDER BY MIRNA_COUNT DESC
  ",
  query_mirna_disease="
    SELECT DISTINCT
    A.MIRNA_PREMATURE AS MIRNA_NAME,
    ANY_VALUE(A.MIRBASE_ACC) AS MIRBASE_ACC,
    C.DISEASE AS DISEASE,
    COUNT(DISTINCT B.PUBMED_ID) AS PUBMED_COUNT,
    GROUP_CONCAT(DISTINCT  B.PUBMED_ID) AS PUBMED_ID
  FROM (select MIRNA_ID,MIRNA_PREMATURE,MIRBASE_ACC from MIRNAS A where 1=1 {mirna_clause}) A
  JOIN MIRNAS_DISEASES_ARTICLES B ON A.MIRNA_ID = B.MIRNA_ID
  JOIN (select DISEASE_ID,DISEASE from DISEASES C WHERE 1=1 {disease_clause}) C ON B.DISEASE_ID = C.DISEASE_ID
  WHERE 1 = 1 {sample_type_clause} {human_clause}
  GROUP BY A.MIRNA_PREMATURE, C.DISEASE
  HAVING COUNT(DISTINCT B.PUBMED_ID) >= {min_assoc_val}
  ORDER BY PUBMED_COUNT DESC
  "
)

query_mg <- list(
  query_by_gene_general="
    SELECT DISTINCT
      B.MIRNA_NAME AS MIRNA_NAME,
      B.MIRBASE_MATURE_ACC AS MIRBASE_MATURE_ACC,
      C.GENE_NAME AS GENE_NAME,
      A.PUBMED_ID AS PUBMED_ID,
      D.PROOF_NAME AS PROOF_NAME
    FROM MIRNAS_GENES_ARTICLES A
    JOIN MIRNAS B ON A.MIRNA_ID = B.MIRNA_ID
    JOIN GENES C ON A.GENE_ID = C.GENE_ID
    JOIN PROOF_MIRNA_GENE D ON A.PROOF_ID = D.PROOF_ID
    {min_assoc_clause}
    WHERE 1=1 {gene_clause_filter} {mirna_clause_filter} {proof_clause} {human_clause}
  ",

  query_by_gene_grouped="
    SELECT DISTINCT
      C.GENE_NAME AS GENE_NAME,
      COUNT(DISTINCT B.MIRNA_NAME) AS MIRNA_COUNT,
      GROUP_CONCAT(DISTINCT B.MIRNA_NAME) AS MIRNA_NAME,
      COUNT(DISTINCT A.PUBMED_ID) AS PUBMED_COUNT,
      GROUP_CONCAT(DISTINCT A.PUBMED_ID) AS PUBMED_ID,
      ANY_VALUE(D.PROOF_NAME) AS PROOF_NAME
    FROM MIRNAS_GENES_ARTICLES A
    JOIN MIRNAS B ON A.MIRNA_ID = B.MIRNA_ID
    JOIN GENES C ON A.GENE_ID = C.GENE_ID
    JOIN PROOF_MIRNA_GENE D ON A.PROOF_ID = D.PROOF_ID
    {min_assoc_clause}
    WHERE 1=1 {gene_clause_filter} {mirna_clause_filter} {proof_clause} {human_clause}
    GROUP BY C.GENE_NAME
    ORDER BY MIRNA_COUNT DESC
  ",

  query_by_gene_plot="
    SELECT DISTINCT
      B.MIRNA_NAME AS MIRNA_NAME,
      C.GENE_NAME AS GENE_NAME,
      COUNT(DISTINCT A.PUBMED_ID) AS PUBMED_COUNT
    FROM MIRNAS_GENES_ARTICLES A
    JOIN MIRNAS B ON A.MIRNA_ID = B.MIRNA_ID
    JOIN GENES C ON A.GENE_ID = C.GENE_ID
    JOIN PROOF_MIRNA_GENE D ON A.PROOF_ID = D.PROOF_ID
    {min_assoc_clause}
    WHERE 1=1 {gene_clause_filter} {mirna_clause_filter} {proof_clause} {human_clause}
    GROUP BY B.MIRNA_NAME, C.GENE_NAME
    ORDER BY PUBMED_COUNT DESC
  ",

  query_by_gene_grouped_mirna="
    SELECT DISTINCT
      B.MIRNA_NAME AS MIRNA_NAME,
      ANY_VALUE(B.MIRBASE_MATURE_ACC) AS MIRBASE_MATURE_ACC,
      COUNT(DISTINCT C.GENE_NAME) AS GENE_COUNT,
      GROUP_CONCAT(DISTINCT C.GENE_NAME) AS GENE_NAME,
      COUNT(DISTINCT A.PUBMED_ID) AS PUBMED_COUNT,
      GROUP_CONCAT(DISTINCT A.PUBMED_ID) AS PUBMED_ID,
      GROUP_CONCAT(DISTINCT D.PROOF_NAME) AS PROOF_NAME
    FROM MIRNAS_GENES_ARTICLES A
    JOIN MIRNAS B ON A.MIRNA_ID = B.MIRNA_ID
    JOIN GENES C ON A.GENE_ID = C.GENE_ID
    JOIN PROOF_MIRNA_GENE D ON A.PROOF_ID = D.PROOF_ID
    {min_assoc_clause}
    WHERE 1=1 {gene_clause_filter} {mirna_clause_filter} {proof_clause} {human_clause}
    GROUP BY B.MIRNA_NAME
    ORDER BY GENE_COUNT DESC
  ",

  query_by_gene_grouped_mirna_gene="
    SELECT DISTINCT
      B.MIRNA_NAME AS MIRNA_NAME,
      ANY_VALUE(B.MIRBASE_MATURE_ACC) AS MIRBASE_MATURE_ACC,
      C.GENE_NAME AS GENE_NAME,
      COUNT(DISTINCT A.PUBMED_ID) AS PUBMED_COUNT,
      GROUP_CONCAT(DISTINCT A.PUBMED_ID) AS PUBMED_ID
    FROM MIRNAS_GENES_ARTICLES A
    JOIN MIRNAS B ON A.MIRNA_ID = B.MIRNA_ID
    JOIN GENES C ON A.GENE_ID = C.GENE_ID
    JOIN PROOF_MIRNA_GENE D ON A.PROOF_ID = D.PROOF_ID
    {min_assoc_clause}
    WHERE 1=1 {gene_clause_filter} {mirna_clause_filter} {proof_clause} {human_clause}
    GROUP BY B.MIRNA_NAME, C.GENE_NAME
    ORDER BY PUBMED_COUNT DESC
  "
)

query_metabolism =list(
  query_db_general="SELECT REACTION_ID, NAME, FORMULA, GPR, SUBSYSTEM, COMPARTMENT FROM REACTIONS ORDER BY NAME;"
)
################################################################################

query_mr=list(
  query_general="SELECT DISTINCT
          B.MIRNA_NAME,
          C.GENE_NAME,
          E.HUMAN_ID,
          E.NAME AS NAME,
          E.SUBSYSTEM,
          E.FORMULA,
          E.GPR,
          E.COMPARTMENT
  FROM (SELECT DISTINCT MIRNA_ID,GENE_ID FROM MIRNAS_GENES_ARTICLES A
        JOIN PROOF_MIRNA_GENE F ON A.PROOF_ID=F.PROOF_ID
        {min_assoc_clause}
        WHERE 1=1 {mirna_clause_filter} {genes_clause_filter} {proof_clause} {human_clause}) A
  JOIN MIRNAS B
  ON A.MIRNA_ID=B.MIRNA_ID
  JOIN (SELECT GENE_ID,GENE_NAME FROM GENES) C
  ON A.GENE_ID=C.GENE_ID
  JOIN REACTIONS_GENES_v D
  ON C.GENE_ID=D.GENE_ID
  JOIN REACTIONS_v E
  ON D.REACTION_ID=E.REACTION_ID",

  # ── Query per ricerca per reazione ─────────────────────────────────────────
  query_by_reaction_general="SELECT DISTINCT
      B.MIRNA_NAME,
      C.GENE_NAME,
      E.HUMAN_ID,
      E.NAME AS NAME,
      E.SUBSYSTEM,
      E.FORMULA,
      E.GPR,
      E.COMPARTMENT
    FROM REACTIONS_v E
    JOIN REACTIONS_GENES_v D ON E.REACTION_ID = D.REACTION_ID
    JOIN GENES C ON D.GENE_ID = C.GENE_ID
    JOIN MIRNAS_GENES_ARTICLES A ON C.GENE_ID = A.GENE_ID
    JOIN MIRNAS B ON A.MIRNA_ID = B.MIRNA_ID
    JOIN PROOF_MIRNA_GENE F ON A.PROOF_ID = F.PROOF_ID
    {min_assoc_clause}
    WHERE 1=1 {reaction_clause_filter} {mirna_clause_filter} {proof_clause} {human_clause}
  "
)


scripts=HTML("
  // Gestione click su link di tipo 'mirna'
  $(document).on('click', '.mirna-link', function(e) {
    e.preventDefault();
    var mirna = $(this).data('mirna');
    Shiny.setInputValue('clicked_mirna', mirna, {priority: 'event'});
  });

  // Gestione click su link di tipo 'disease'
  $(document).on('click', '.disease-link', function(e) {
    e.preventDefault();
    var disease = $(this).data('disease');
    Shiny.setInputValue('clicked_disease', disease, {priority: 'event'});
  });

  // Gestione click su link di tipo 'pubmed'
  $(document).on('click', '.pubmed-link', function(e) {
    e.preventDefault();
    var pubmed = $(this).data('pubmed');
    Shiny.setInputValue('clicked_pubmed', pubmed, {priority: 'event'});
  });
  
  // Gestione click su link di tipo 'pubmed'
  $(document).on('click', '.gene-link', function(e) {
    e.preventDefault();
    var gene = $(this).data('gene');
    Shiny.setInputValue('clicked_gene', gene, {priority: 'event'});
  });
  

")
  

documentation_miRNA_disease_table_upload = HTML("
  <b>Upload your own miRNA–disease associations.</b> The original database is never overwritten —
  your data replaces it only for the current session and can be restored at any time.<br><br>
  <b>Required columns:</b><br>
  &bull; <code>MIRNA_NAME</code> — premature miRNA name (e.g. <i>hsa-miR-21</i>)<br>
  &bull; <code>DISEASE</code> — disease name matching the Disease Ontology (e.g. <i>lung cancer</i>)<br>
  &bull; <code>PUBMED_ID</code> — integer PubMed ID<br><br>
  <b>Optional columns:</b><br>
  &bull; <code>SAMPLE_TYPE</code> — one of: circulating, exosome, tissue, other<br><br>
  Accepted formats: <code>.csv</code>, <code>.tsv</code>, <code>.xlsx</code>
")
documentation_metabolism_table_upload = HTML("
  <b>Upload a custom metabolic model.</b> The original database is never overwritten —
  your data replaces it only for the current session and can be restored at any time.<br><br>
  <b>Required columns:</b><br>
  &bull; <code>NAME</code> — reaction name<br>
  &bull; <code>FORMULA</code> — reaction formula (e.g. <i>A + B =&gt; C</i>)<br>
  &bull; <code>GPR</code> — gene–protein–reaction rule (e.g. <i>GENE1 or GENE2</i>)<br>
  &bull; <code>SUBSYSTEM</code> — metabolic subsystem / pathway<br>
  &bull; <code>COMPARTMENT</code> — cellular compartment<br><br>
  <b>Optional columns:</b><br>
  &bull; <code>GENES</code> — comma-separated gene symbols; required to update gene–reaction associations<br>
  &bull; <code>HUMAN_ID</code> — original reaction identifier (e.g. <i>MAR03905</i>)<br><br>
  Accepted formats: <code>.csv</code>, <code>.tsv</code>, <code>.xlsx</code>
")

# =============================================================================
# Subsystem → biological category mapping (Human-GEM, 138 subsystems → 10 categories)
# NA_character_ entries are excluded from category view (transport, drug, misc.)
# =============================================================================
SUBSYSTEM_CATEGORY_MAP <- c(
  # ── Fatty acid metabolism ─────────────────────────────────────────────────
  "Fatty acid oxidation"                                              = "Fatty acid metabolism",
  "Fatty acid activation (cytosolic)"                                 = "Fatty acid metabolism",
  "Fatty acid activation (endoplasmic reticular)"                     = "Fatty acid metabolism",
  "Carnitine shuttle (cytosolic)"                                     = "Fatty acid metabolism",
  "Carnitine shuttle (mitochondrial)"                                 = "Fatty acid metabolism",
  "Carnitine shuttle (endoplasmic reticular)"                         = "Fatty acid metabolism",
  "Carnitine shuttle (peroxisomal)"                                   = "Fatty acid metabolism",
  "Fatty acid biosynthesis"                                           = "Fatty acid metabolism",
  "Fatty acid biosynthesis (even-chain)"                              = "Fatty acid metabolism",
  "Fatty acid biosynthesis (unsaturated)"                             = "Fatty acid metabolism",
  "Fatty acid biosynthesis (odd-chain)"                               = "Fatty acid metabolism",
  "Fatty acid metabolism"                                             = "Fatty acid metabolism",
  "Fatty acid elongation (even-chain)"                                = "Fatty acid metabolism",
  "Fatty acid elongation (odd-chain)"                                 = "Fatty acid metabolism",
  "Fatty acid desaturation (even-chain)"                              = "Fatty acid metabolism",
  "Fatty acid desaturation (odd-chain)"                               = "Fatty acid metabolism",
  "Acyl-CoA hydrolysis"                                               = "Fatty acid metabolism",
  "Acylglycerides metabolism"                                         = "Fatty acid metabolism",
  "Fatty acid degradation"                                            = "Fatty acid metabolism",

  # ── Beta-oxidation ────────────────────────────────────────────────────────
  "Beta oxidation of even-chain fatty acids (mitochondrial)"          = "Beta-oxidation",
  "Beta oxidation of even-chain fatty acids (peroxisomal)"            = "Beta-oxidation",
  "Beta oxidation of odd-chain fatty acids (mitochondrial)"           = "Beta-oxidation",
  "Beta oxidation of unsaturated fatty acids (n-9) (mitochondrial)"  = "Beta-oxidation",
  "Beta oxidation of unsaturated fatty acids (n-9) (peroxisomal)"    = "Beta-oxidation",
  "Beta oxidation of unsaturated fatty acids (n-7) (mitochondrial)"  = "Beta-oxidation",
  "Beta oxidation of unsaturated fatty acids (n-7) (peroxisomal)"    = "Beta-oxidation",
  "Beta oxidation of poly-unsaturated fatty acids (mitochondrial)"   = "Beta-oxidation",
  "Beta oxidation of di-unsaturated fatty acids (n-6) (mitochondrial)" = "Beta-oxidation",
  "Beta oxidation of di-unsaturated fatty acids (n-6) (peroxisomal)" = "Beta-oxidation",
  "Beta oxidation of branched-chain fatty acids (mitochondrial)"      = "Beta-oxidation",
  "Beta oxidation of phytanic acid (peroxisomal)"                     = "Beta-oxidation",
  "Omega-3 fatty acid metabolism"                                     = "Beta-oxidation",
  "Omega-6 fatty acid metabolism"                                     = "Beta-oxidation",
  "Linoleate metabolism"                                              = "Beta-oxidation",

  # ── Eicosanoid / Inflammatory lipids ─────────────────────────────────────
  "Leukotriene metabolism"                                            = "Eicosanoid / Inflammatory lipids",
  "Arachidonic acid metabolism"                                       = "Eicosanoid / Inflammatory lipids",
  "Prostaglandin biosynthesis"                                        = "Eicosanoid / Inflammatory lipids",
  "Eicosanoid metabolism"                                             = "Eicosanoid / Inflammatory lipids",

  # ── Sterol metabolism ─────────────────────────────────────────────────────
  "Bile acid biosynthesis"                                            = "Sterol metabolism",
  "Formation and hydrolysis of cholesterol esters"                    = "Sterol metabolism",
  "Bile acid recycling"                                               = "Sterol metabolism",
  "Cholesterol metabolism"                                            = "Sterol metabolism",
  "Cholesterol biosynthesis 1 (Bloch pathway)"                        = "Sterol metabolism",
  "Cholesterol biosynthesis 2"                                        = "Sterol metabolism",
  "Cholesterol biosynthesis 3 (Kandustch-Russell pathway)"            = "Sterol metabolism",
  "Steroid metabolism"                                                = "Sterol metabolism",
  "Androgen metabolism"                                               = "Sterol metabolism",
  "Estrogen metabolism"                                               = "Sterol metabolism",
  "Glucocorticoid biosynthesis"                                       = "Sterol metabolism",
  "Vitamin D metabolism"                                              = "Sterol metabolism",
  "Vitamin E metabolism"                                              = "Sterol metabolism",
  "Retinol metabolism"                                                = "Sterol metabolism",
  "Vitamin A metabolism"                                              = "Sterol metabolism",

  # ── Complex lipids ────────────────────────────────────────────────────────
  "Sphingolipid metabolism"                                           = "Complex lipids",
  "Glycerophospholipid metabolism"                                    = "Complex lipids",
  "Glycerolipid metabolism"                                           = "Complex lipids",
  "Glycosphingolipid biosynthesis-lacto and neolacto series"          = "Complex lipids",
  "Glycosphingolipid biosynthesis-ganglio series"                     = "Complex lipids",
  "Glycosphingolipid biosynthesis-globo series"                       = "Complex lipids",
  "Glycosphingolipid metabolism"                                      = "Complex lipids",
  "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis"            = "Complex lipids",
  "Inositol phosphate metabolism"                                     = "Complex lipids",
  "Phosphatidylinositol phosphate metabolism"                         = "Complex lipids",
  "Ether lipid metabolism"                                            = "Complex lipids",

  # ── Amino acid metabolism ─────────────────────────────────────────────────
  "Phenylalanine, tyrosine and tryptophan biosynthesis"               = "Amino acid metabolism",
  "Arginine and proline metabolism"                                   = "Amino acid metabolism",
  "Glycine, serine and threonine metabolism"                          = "Amino acid metabolism",
  "Valine, leucine, and isoleucine metabolism"                        = "Amino acid metabolism",
  "Tyrosine metabolism"                                               = "Amino acid metabolism",
  "Alanine, aspartate and glutamate metabolism"                       = "Amino acid metabolism",
  "Cysteine and methionine metabolism"                                = "Amino acid metabolism",
  "Lysine metabolism"                                                 = "Amino acid metabolism",
  "Histidine metabolism"                                              = "Amino acid metabolism",
  "Tryptophan metabolism"                                             = "Amino acid metabolism",
  "Phenylalanine metabolism"                                          = "Amino acid metabolism",
  "Serotonin and melatonin biosynthesis"                              = "Amino acid metabolism",
  "Aminoacyl-tRNA biosynthesis"                                       = "Amino acid metabolism",
  "Beta-alanine metabolism"                                           = "Amino acid metabolism",
  "Metabolism of other amino acids"                                   = "Amino acid metabolism",

  # ── Nucleotide metabolism ─────────────────────────────────────────────────
  "Nucleotide metabolism"                                             = "Nucleotide metabolism",
  "Purine metabolism"                                                 = "Nucleotide metabolism",
  "Pyrimidine metabolism"                                             = "Nucleotide metabolism",
  "Amino sugar and nucleotide sugar metabolism"                       = "Nucleotide metabolism",

  # ── Carbohydrate / Energy ─────────────────────────────────────────────────
  "Glycolysis / Gluconeogenesis"                                      = "Carbohydrate / Energy",
  "Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism"  = "Carbohydrate / Energy",
  "Pyruvate metabolism"                                               = "Carbohydrate / Energy",
  "Pentose phosphate pathway"                                         = "Carbohydrate / Energy",
  "Starch and sucrose metabolism"                                     = "Carbohydrate / Energy",
  "Fructose and mannose metabolism"                                   = "Carbohydrate / Energy",
  "Galactose metabolism"                                              = "Carbohydrate / Energy",
  "Oxidative phosphorylation"                                         = "Carbohydrate / Energy",
  "Butanoate metabolism"                                              = "Carbohydrate / Energy",
  "Propanoate metabolism"                                             = "Carbohydrate / Energy",
  "Pentose and glucuronate interconversions"                          = "Carbohydrate / Energy",
  "C5-branched dibasic acid metabolism"                               = "Carbohydrate / Energy",

  # ── Glycan / ECM ──────────────────────────────────────────────────────────
  "Keratan sulfate degradation"                                       = "Glycan / ECM",
  "Keratan sulfate biosynthesis"                                      = "Glycan / ECM",
  "N-glycan metabolism"                                               = "Glycan / ECM",
  "Chondroitin / heparan sulfate biosynthesis"                        = "Glycan / ECM",
  "Chondroitin sulfate degradation"                                   = "Glycan / ECM",
  "Heparan sulfate degradation"                                       = "Glycan / ECM",
  "O-glycan metabolism"                                               = "Glycan / ECM",
  "Blood group biosynthesis"                                          = "Glycan / ECM",

  # ── Cofactors / Vitamins ──────────────────────────────────────────────────
  "Folate metabolism"                                                 = "Cofactors / Vitamins",
  "Biopterin metabolism"                                              = "Cofactors / Vitamins",
  "Nicotinate and nicotinamide metabolism"                            = "Cofactors / Vitamins",
  "Vitamin B6 metabolism"                                             = "Cofactors / Vitamins",
  "Biotin metabolism"                                                 = "Cofactors / Vitamins",
  "Thiamine metabolism"                                               = "Cofactors / Vitamins",
  "Riboflavin metabolism"                                             = "Cofactors / Vitamins",
  "Vitamin B2 metabolism"                                             = "Cofactors / Vitamins",
  "Lipoic acid metabolism"                                            = "Cofactors / Vitamins",
  "Vitamin B12 metabolism"                                            = "Cofactors / Vitamins",
  "Pantothenate and CoA metabolism"                                   = "Cofactors / Vitamins",
  "Terpenoid backbone biosynthesis"                                   = "Cofactors / Vitamins",
  "Ubiquinone and other terpenoid-quinone biosynthesis"               = "Cofactors / Vitamins",
  "Ubiquinone synthesis"                                              = "Cofactors / Vitamins",
  "Porphyrin metabolism"                                              = "Cofactors / Vitamins",
  "Heme degradation"                                                  = "Cofactors / Vitamins",
  "Glutathione metabolism"                                            = "Cofactors / Vitamins",
  "ROS detoxification"                                                = "Cofactors / Vitamins",
  "Sulfur metabolism"                                                 = "Cofactors / Vitamins",
  "Ascorbate and aldarate metabolism"                                 = "Cofactors / Vitamins",
  "Vitamin C metabolism"                                              = "Cofactors / Vitamins",
  "Glyoxylate and dicarboxylate metabolism"                           = "Cofactors / Vitamins",

  # ── Protein metabolism ────────────────────────────────────────────────────
  "Protein degradation"                                               = "Protein metabolism",
  "Protein assembly"                                                  = "Protein metabolism",
  "Protein modification"                                              = "Protein metabolism",
  "Peptide metabolism"                                                = "Protein metabolism",

  # ── Excluded (transport, drug, misc — hidden in category view) ───────────
  "Transport reactions"                                               = "Transport reactions",
  "Drug metabolism"                                                   = NA_character_,
  "Xenobiotics metabolism"                                            = NA_character_,
  "Isolated"                                                          = NA_character_,
  "Miscellaneous"                                                     = NA_character_,
  "Pool reactions"                                                    = NA_character_,
  "Insect hormone biosynthesis"                                       = NA_character_,
  "Toluene degradation"                                               = NA_character_,
  "octane oxidation"                                                  = NA_character_
)
