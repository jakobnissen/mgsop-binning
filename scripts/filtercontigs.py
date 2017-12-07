# Renames contigs to unique names and removes contigs below minimum size.

# Assigning this value to local namespace avoids dotted lookup
min_contig_length = snakemake.params.min_contig_length

def iterfasta(filehandle):
    "A generator which yield header, seq tuples from an open fasta file."
    
    buffer = list()
    
    header = next(filehandle).strip()
    if not header.startswith('>'):
        raise ValueError('First line is not a header')
    header = header[1:]
    
    for line in map(str.rstrip, filehandle):
        if line.startswith('>'): 
            yield header, ''.join(buffer)
            buffer.clear()
            header = line[1:]
            
        else:
            buffer.append(line)
            
    yield header, ''.join(buffer)

with open(snakemake.output[0], 'w') as outfile:
    for filename in snakemake.input:
        samplename = filename.rpartition('/')[2].partition('.')[0]

        with open(filename) as infile:
            for header, seq in iterfasta(infile):
                if len(seq) >= min_contig_length:
                    print('>{}_{}'.format(samplename, header), file=outfile)
                    print(seq, file=outfile)
