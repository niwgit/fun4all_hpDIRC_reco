# fun4all_hpDIRC_reco

git clone https://github.com/niwgit/fun4all_hpDIRC_reco

Run geometric reconstruction macro:
root -l 'geo_reco_f4a.cpp("filename.root","lut_avr.root",2)'

Run time imaging reconstruction macros:
1) To create pdfs
root -l 'createPdf_f4a.cpp("filename.root")'

2) For reconstruction
root -l 'recoPdf_f4a.cpp("filename.root","filename.pdf.root")'

where "filename.root" is the input root file. 