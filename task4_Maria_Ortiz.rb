require 'bio'
require 'stringio'
require 'io/console'

#"This property may, in part, be due to their putative regulatory
#role as the dominant reading frame shifts and length of the overlaps suggests.
#Overlapping genes in Prochlorococcus spp.,which are closely phylogenetically
#related yet have significantly different nucleotide sequences for orthologs, further support this role"

#BIBLIOGRAPHY: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC525685/

#I took 10^-6 e - values and more of 50% of coverage (overlapping) 
#https://academic.oup.com/bioinformatics/article/24/3/319/252715

#Create a new file
def new_file(filename)
  return (File.open(filename, "w"))
end


#Fasta files and type of the files for specie1 and specie2
specie1_file =  "pep.fa"
specie2_file = "TAIR10_seq_20110103_representative_gene_model_updated"


#New file output 
output= new_file("output_orthologues_overlap_50_evalue_10-6.txt")
output.puts "Assignment 4 - Maria Ortiz Rodriguez"
output.puts "ORTHOLOGUES FOUND in files #{specie1_file} and #{specie2_file}\n"


# Creation of the folders
system("rm -r Databases")
system("mkdir Databases")


# Creation of databases
db_specie1 = (specie1_file.to_s + '_db')
db_specie2 = (specie2_file.to_s + '_db')


# We create the databases
system("makeblastdb -in '#{specie1_file}' -dbtype prot -out ./Databases/#{db_specie1.to_s}") 
system("makeblastdb -in '#{specie2_file}' -dbtype nucl -out ./Databases/#{db_specie2.to_s}")


# We create factories depending on the types of the files
factory_specie1 = Bio::Blast.local('blastx', "./Databases/#{db_specie1.to_s}")#, -e 10**-6)
factory_specie2 = Bio::Blast.local('tblastn', "./Databases/#{db_specie2.to_s}")# -e 10**-10)


# Open fasta file
fasta_specie1 = Bio::FastaFormat.open(specie1_file)
fasta_specie2 = Bio::FastaFormat.open(specie2_file)

#Create a hash with sequences and id of the second species
hash_specie2 = Hash.new
fasta_specie2.each do |feature|
  hash_specie2[feature.entry_id] = (feature.seq).to_s
end

#Blast
fasta_specie1.each do |seq1|
  
  specie1_id = seq1.entry_id.to_s #first specie id
  
  query1 = factory_specie2.query(seq1)
  
  unless query1.hits[0].nil?
    
    unless query1.hits[0].evalue.nil?
      
      overlap1 = query1.hits[0].overlap.to_i
      evalue1 = query1.hits[0].evalue.to_i

      
      if (evalue1<= 10**(-6)) && (overlap1 >= 50) #percentage to be homologous
        specie2_id = query1.hits[0].definition.match(/^\w+.\w+/).to_s
        
        ################################################################################
        
        if hash_specie2.has_key?(specie2_id)
        
        result1 = factory_specie1.query(">#{specie2_id}\n#{hash_specie2[specie2_id]}") 
        
          unless result1.hits[0].nil?
          
            unless result1.hits[0].evalue.nil?

                overlap2 = (result1.hits[0].overlap.to_i)
                evalue2 = (result1.hits[0].evalue.to_i)
                
              if (evalue2 <= 10**(-6)) && (overlap2 >=50)
#percentage to be orthologs - this is not enough to filter for orthologs but watching the scores its enough to filter the most significant in this case 

                match = result1.hits[0].definition.match(/^\w+.\w+/).to_s
              
                if specie1_id == match
                output.puts "#{specie1_id} -- #{specie2_id}"
                end
              
              end
            
            end#segundo unless
        
          end #prmer unless
          
        end #if hash has key
      
      end
      
    end
    
  end

end



output.puts "\n\nBonus exercise: \nWe can't decide if those species are orthologs or not just with this analysis. To improve the knowledge between this species we can also make a phylogenetic tree"


output.close
system("rm -r Databases") 

