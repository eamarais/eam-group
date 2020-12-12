##############################################################################################
##############################################################################################
# This script provides solutions to problems in handling pdfs which are freqently met in research
# merge PDFs and extract certain pages

# Motivations of using Python:
# free and safe
# fast for handling massive files
# flexible and powerful (edit your own codes for your specific tasks)

# You can use functions from this script: from handle_pdfs_for_research import "xxx" 
# but make sure the "PyPDF2" package is installed in advance 
# Or you can use the script interactively and develop what you need
##############################################################################################
# set up working directory
import os 
os.chdir('C:\\Users\\geosc\\Python_workshop')

# get the list of pdf files
pdfs_list = ['test1.pdf','test2.pdf','test3.pdf']
##############################################################################################
##############################################################################################
# fundamental solutions to most PDF problems (merge PDFs and extract certain pages)

# create a function to merge a list of pdfs in order (merge process = "append" + "save out")
# usage: merge_pdfs(pdfs_list,"combine_pdf_test1")

def merge_pdfs(pdfs_list,output_filename):
    '''Input the filenames of a list of pdfs, combine them into a single pdf with the specified name, 
       keeping the default input order.'''  
    # initialize the tool for merging pdf ("merger" here)
    from PyPDF2 import PdfFileMerger
    merger = PdfFileMerger()
    # append each input pdf together
    for pdf in pdfs_list:
        merger.append(pdf)
    # save out the combined file
    merger.write(str(output_filename)+'.pdf')
    merger.close()

# create a function to split a multi-pages pdf and save out each selected page as a separate pdf
# usage: split_pdfs('test1.pdf',"splite_pdf_test1",[1,2,5])

def split_pdfs(input_pdf,output_filename,target_pages=None):
    '''Input a single pdf filename, save out each page as a separate pdf by default.
       Or you can provide selected page numbers in a list (just use the natural page number in the file).
       Outpout filenames are "Your filename + page number"'''
    # import the tools needed
    from PyPDF2 import PdfFileWriter, PdfFileReader
    # read the input pdf
    inputpdf = PdfFileReader(open(input_pdf, "rb"))
    # save out each page if no arguments are provided
    # or save out the requested target pages
    # when providing the target pages, just use the natural page number in the file
    # they are converted to Python indice in the codes
    if target_pages is None:
        output_pages = range(inputpdf.numPages)
    else: 
        output_pages = [x - 1 for x in target_pages]
    # loop over all pages or requested ones only, then save out each page
    for i in output_pages:
        output = PdfFileWriter()
        output.addPage(inputpdf.getPage(i))
        with open(str(output_filename)+'_page_%s.pdf' % int(i+1), "wb") as outputStream:
             output.write(outputStream)
##############################################################################################
# Example application: merge all pdf files within a directory

# list all files within this directory
all_files = os.listdir()

# select your pdfs
# there are many ways to select the file, the key is to play around the "string"
# first get the string for each file name, then use "str.startswith", "str.endswith"
# or use if statements like: if "xxx" in str

pdf_files = []

# loop over all files
for file in all_files:
    # get the strings for each file
    filename = os.fsdecode(file) 
    #  only add the file to your list if it is a pdf
    if filename.endswith('.pdf') == True: 
        pdf_files.append(filename)

# sort the images to decide the order 
pdf_files.sort() 

# check the files and the order
for file in pdf_files:
    print(file,sep='\n')
    
# merge together
merge_pdfs(pdf_files,"final_combination")
##############################################################################################
# Other PDF merging problems are not targted for now, as they are not very frequent
# And anyway, using "merge_pdfs" and "split_pdfs" together can already solve a range of problems
# Here I only list some other potentials from this package or Python in general

# only combine the first 3 pages from each file
merger = PdfFileMerger()

for pdf in pdfs:
    merger.append(pdf,pages=(0,3)) 
    
merger.write("combined_pdf_test3.pdf")
merger.close()
  
# only combine pages 1,3,5 from each input pdf file
merger = PdfFileMerger()

for pdf in pdfs:
    merger.append(pdf, pages=(0,6,2)) # (start,end,sep)
    
merger.write("combined_pdf_test4.pdf")
merger.close()

# merge two pages into one
# https://stackoverflow.com/questions/22795091/how-to-append-pdf-pages-using-pypdf2

# End
##############################################################################################
##############################################################################################
