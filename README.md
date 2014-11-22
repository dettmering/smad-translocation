smad-translocation
==================

CellProfiler and R analysis pipeline for determining Smad protein translocation into the nucleus. This pipeline detects cytoplasm (i.e. stained with Phalloidin-PI) and nucleus (stained with DAPI) of each cell in the image. A third channel is used to detect Smad translocation by comparing the median fluorescence intensity of cytoplasm and nucleus per cell.

You might want to compensate for uneven background illumination, this is not included in the pipeline.
