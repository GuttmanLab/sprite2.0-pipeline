import sys
import operator

fileName=sys.argv[1]
maxSize=sys.argv[2]
gene=sys.argv[3]
qChrom=sys.argv[4]
qStart=sys.argv[5]
qEnd=sys.argv[6]
weighted=sys.argv[7]
plusMin=sys.argv[8]

linesClusters = open(fileName).read().splitlines()#[1:100]

#make a dictionary 
bins={}
chromosomes = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19','chrX']
res=100000
for line in linesClusters:
	#print line
	barcode= line.split()[0]
	#print(barcode)
	clusters=line.split()[1:]
	#print(clusters)
	
	if len(clusters) > int(maxSize):
		continue
	xist = False
	for molecule in clusters:
		tag=molecule.split("_")[0].split("(")[0]
		strand=molecule.split("_")[0].split("(")[1].split(")")[0]
		#print(tag,strand)
		chrom="chr"+molecule.split("_")[1].split(":")[0]
		#print(chrom)
		pos=int(molecule.split(":")[1])
		start=int(pos/res)
		#print(pos,start)

		if tag == "RPM" and chrom == qChrom and strand == plusMin:
			if int(qStart) < pos < int(qEnd):
				xist = True
				#print(molecule,plusMin)
	
	if weighted == "n_over_two":
		weight = 2/len(clusters)
	else:
		weight = 1

	if xist == True:
		for molecule in clusters:
			tag=molecule.split("_")[0].split("(")[0]
			if tag == "RPM":
				continue
			#print(molecule)
			chrom="chr"+molecule.split("_")[1].split(":")[0]
			pos=int(molecule.split(":")[1])
			start=int(pos/res)
			#print(chrom,pos,start)
			if chrom not in bins:
				bins[chrom] = {}
			if start not in bins[chrom]:
				bins[chrom][start]=0
			bins[chrom][start] = bins[chrom][start]+weight

#print(bins)
fout=open(gene+"_"+str(maxSize)+"_res-"+str(res)+"_weight-"+weighted+".bedgraph",'w')
for i in range(0, len(chromosomes)):
	chr = chromosomes[i]
	
	try:
		a=bins[chr]
	except KeyError:
		continue
	for a,b in sorted(bins[chr].items(),key = operator.itemgetter(0), reverse = False):
		if a in bins[chr]:
			fout.write(str(chr)+"\t"+str(int(a*res))+"\t"+str(int(a*res+res-1))+"\t"+str(b)+"\n")
fout.close()

#for chromo in chromosomes:	
#	try:
#		for a,b in sorted(bins[chromo].items(),key = operator.itemgetter(0), reverse = False):
#			fout.write(str(chromo)+"\t"+str(a*res)+"\t"+str(a+res)+ "\t"+str(b)+"\n")
#	except KeyError:
#		continue
#fout.close()
