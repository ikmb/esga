#!/usr/bin/env ruby
# A script to download the latest version of an EnsEMBL genome,
# including GTF files.
# = USAGE
# ./get_ensembl_genome.rb gorilla_gorilla 67

require 'fileutils'
require 'net/ftp'
require 'optparse'
require 'ostruct'

def exec(command)

	warn "Running: #{command}"
	system(command)
end

options = OpenStruct.new()
opts = OptionParser.new()
opts.on("-v","--version", "=VERSION","EnsEMBL version") {|argument| options.version = argument }
opts.on("-s","--species", "=SPECIES","Species name") {|argument| options.species = argument }
opts.on("-o","--outfile", "=OUTFILE","Output") {|argument| options.outfile = argument }
opts.on("-h","--help","Display the usage information") {
  puts opts
  exit
}

opts.parse! 

species = options.species or fail "Must provide species in snake_case (-s)"
version = options.version or fail "Must provide release version (-v) to use"

species.downcase!

#raise "This does not look like a valid species name (must be lower case scientific name with underscore)" unless species.match(/[a-z]*_[a-z]*/)

# Rsync URL for EnsEMBL 
RSYNC_PATH = "rsync://ftp.ensembl.org/ensembl/pub/release-#{version}"
# FTP Path for EnsEMBL
WGET_PATH = "ftp://ftp.ensembl.org/pub/release-#{version}"

ftp=Net::FTP.new
ftp.connect("ftp.ensembl.org")
ftp.passive = true
ftp.login

# check presence of directories first

species.split(":").each do |s|

	begin
		puts "Downloading genome for #{s}, release #{version}"
		ftp.chdir("/pub/release-#{version}/fasta/#{s}/dna/")

		file = ftp.nlst.find{|f| f.include?("dna_rm.primary_assembly.fa.gz") }
	
		file = ftp.nlst.find{|f| f.include?("dna_rm.toplevel.fa.gz")} if file.nil? or file.length < 1
	
		warn "Downloading genome sequence: #{file}"
		ftp.get(file,File.basename(file))

		warn "Extracting genomes"
		exec("gunzip -c #{file} | sed 's/ .*//'  > #{s}.fa")
		exec("rm #{file}")
		

	rescue
		raise "Something went wrong. This species name was probably not found on the server (Check spelling!)"
	end

	begin
		warn "Downloading GTF for #{}"
		ftp.chdir("/pub/release-#{version}/gtf/#{s}")
		file = ftp.nlst.find{|f| f.include?("#{version}.gtf.gz")}
		raise "No GTF file was found for this species" if file.nil?
	
		#ftp.get(file,File.basename(file))
		exec("wget #{WGET_PATH}/gtf/#{s}/#{file}")
	
		exec("gunzip -c #{file} > #{s}.gtf")
		exec("rm #{file}")

	rescue
		raise "Something went wrong. This species name was probably not found on the server (Check spelling!)"
	end

end

ftp.close

warn "All files have been retrieved and formatted"



