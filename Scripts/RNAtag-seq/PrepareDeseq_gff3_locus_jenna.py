from sys import argv
script, input_file = argv

new_list = []

gfflist = open(input_file, "r").readlines()
for line1 in gfflist:
    identifier = line1.rstrip().split('\t')[2]
    if identifier in ["gene", "ncRNA", "pseudogene"]:
       attribute = line1.rstrip().split('\t')[8]
       #print attribute
       if "locus_tag" in attribute:
          locus = attribute.split('locus_tag=')[-1]
          reads = line1.rstrip().split('\t')[9]
          #print locus
       else:
          locus = attribute.split(';')[1].split('=')[1] + "_" +attribute.split(';')[0].split('=')[1]
          reads = line1.rstrip().split('\t')[9]
          #locus.split("=")[-1] == "tRNA":
          
          print locus
       

    #print locus
           #print locus_final
       new_list.append([locus, reads])

out_file = open(input_file.replace(".txt", "_filt.txt"), "w")
out_line = ["%s\t%s" % (a[0], a[1]) for a in new_list]
out_file.write('\n'.join(out_line))
out_file.close()
