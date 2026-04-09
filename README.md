# Study Discovery Engine

An AI-powered research tool that automatically finds, filters, and 
evaluates public microbiome studies for inclusion in large-scale 
meta-analyses like the Gut Microbiome Tree of Life project (GMToL).

## What it does

1. **Searches** SRA and ENA programmatically for candidate studies
   using taxonomy-based queries
2. **Filters** studies against hard criteria — paired 16S and WGS 
   data, downloadable FASTQs, valid accession codes
3. **Scores** candidate studies by scientific value — phylogenetic 
   gap coverage, metadata completeness, geographic representation
4. **Presents** findings via an AI agent that summarizes each study,
   assesses its value for GMToL, and asks targeted questions
5. **Learns** from researcher feedback to improve recommendations

## Architecture
SRA/ENA/QIITA --> src/fetcher.py (fetch and paginate study metadata) --> src/filters.py (hard filter checks such as paired 16S + WGS, has FASTQs etc.) --> AI agent (summarizes, assesses, asks questions) --> Researcher feedback (gets added to decision log)
