CC = gcc
CFLAGS = -g -Wall -O2
LDFLAGS = -lhts
TARGET = get_consensus_fasta_from_vcf

all: $(TARGET)

$(TARGET): get_consensus_fasta_from_vcf.c
	$(CC) $(CFLAGS) -o $(TARGET) $< $(LDFLAGS)


#The test VCF contains two records for two individuals with the following data:
#	Record 1:
#	- Sample ind1: 0|0 (homozygous, REF allele)
#	- Sample ind2: 0|1 (heterozygous; random choice expected to pick ALT)
#	Record 2:
#	- Sample ind1: 0|1 (heterozygous; random choice expected to pick ALT)
#	- Sample ind2: 0|0 (homozygous, REF allele)
#With a fixed seed (42), we expect the consensus FASTA output to be:
#	>ind1
#	AC
#	>ind2
#	AC
test: $(TARGET)
	./$(TARGET) 42 test1.vcf test1.fa
	diff test1.fa test1_expected.fa

clean:
	rm -f $(TARGET) *.o