#%%
transcripts = {}
# unique = set()
current_transcript = ''
transcript_order = []
transcript_info = {}
with open('/mnt/data/genomes/GRCm38.p6/GRCm38.p6.annotation.gtf') as f:
	for line in f:
		if len(line.split()) > 13:
			if line.split()[2] == 'transcript':
				chrom = line.split()[0]
				strand = line.split()[6]
				start = int(line.split()[3])
				end = int(line.split()[4])
				current_transcript = line.split()[13].strip(';').strip('"')
				transcript_info[current_transcript] = line.strip()
				transcripts[current_transcript] = [ chrom, strand, [start, end] ]
				transcript_order.append(current_transcript)
			if line.split()[2] == 'exon':
				start = int(line.split()[3])
				end = int(line.split()[4])
				transcripts[current_transcript].append( [start,end] )

#%%
with open('/mnt/data/genomes/GRCm38.p6/GRCm38.p6.exon_intron.gtf', 'w') as out:		
	for transcript in transcript_order:
		data = transcripts[transcript]
		chrom = data[0]
		strand = data[1]
		start = data[2][0]
		end = data[2][1]
		exons = data[3:]
		if strand == '-':
			exons.reverse()
		out.write(transcript_info[transcript] + '\n')
		type = transcript_info[transcript].split()[1]
		rest_tab = '\t'.join(transcript_info[transcript].split()[5:8])
		rest = ' '.join(transcript_info[transcript].split()[8:])
		for i in range(0, len(exons)):# in exons:
			if i == 0:
				out.write(chrom+'\t'+type+'\texon\t'+str(exons[i][0])+'\t'+str(exons[i][1])+'\t'+rest_tab +'\t'+rest + '\n')
			else:
				out.write(chrom+'\t'+type+'\tintron\t'+str(exons[i-1][1]+1)+'\t'+str(exons[i][0]-1)+'\t'+rest_tab +'\t'+rest + '\n')
				out.write(chrom+'\t'+type+'\texon\t'+str(exons[i][0])+'\t'+str(exons[i][1])+'\t'+rest_tab+'\t'+rest + '\n')

#%%
