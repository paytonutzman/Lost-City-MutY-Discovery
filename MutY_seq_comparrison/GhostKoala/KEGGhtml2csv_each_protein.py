# modification of original KEGGhtml2csv.py to report each protein (KEGG protein IDs) instead of each module

import glob
for f in glob.glob("*.htm*"):
    d = {}
    previous = 'blank'
    outfilename = f.replace('.html','.csv').replace('.htm','.csv')
    outfilename = outfilename.replace('.csv','.proteins.csv')
    with open (f, 'rt') as myfile:
        for line in myfile:
            if 'target="_blank">K' in line:
            	protid = line.split('target="_blank">')
            	protid = protid[1].split('</a></dt>')
            	protid = protid[0]
            	d[protid] = ''			# blank entry in the dictionary for now
            	previous = protid
            elif '<dd>' in line:
	            if not previous == 'blank': 	
            		contigs = line.split('<dd>')
            		contigs = contigs[1].split('</dd>')
            		contigs = contigs[0].replace(',','.')	# replace commas with periods to avoid csv formatting failure
            		d[previous] = contigs	# assigns contigs to the previous protein ID
            
        
    with open(outfilename, "w") as outfile:		
        outfile.write('Protein,')
        outfile.write(outfilename.replace('htmldir/','').replace('.csv',''))
        outfile.write('\n')
        for i in d:
            outfile.write(i + ",")
            outfile.write(d[i])
            outfile.write("\n")
            
