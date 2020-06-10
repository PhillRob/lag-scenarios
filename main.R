
infile = "NZ-species-fre-input.csv" #Frequency input data
infile2 = "NZ-species-annual-fre-input.csv"    # Annual number of specimens collected data
outfile = 'NZ_Species_Output.csv'   #Output file

## Read the data

data = read.csv(infile, header=T)

# Prepare the data for the main part of the code
# It is important to have the columns named "Year", Frequency" and "Specimens" to get it work
# If your data already has this, you may omit

names(data)[5] = "Frequency"  # 5 th column of the input contains Frequncy
names(data)[7] = "Specimens"  # 7th column of the input file contains no_specimens
data$Species = as.character(data$Species)

## Read the data where only specimens are there but no frequency count

ydata = read.csv(infile2, header=T)
names(ydata)[2] = "Year"
names(ydata)[3] = "Specimens"

# Load all the functions

source("utilities.R")

### Main Part of the Code ###

## Find out all the species and fit the model for each of them
## We run separately for including zeros and exclusing zeros

out0 = lagfit(data, ydata, zeros=TRUE) # With all zeros
out1 = lagfit(data, ydata, zeros=FALSE) # Without zeros

names(out1)[2:6] = paste(names(out1)[3:7],"_nonzero") #Add Nonzero to the labels of nonzero fits

# Merge all of them in a bigger data frame to write

outdata = merge(out0, out1, by=c("Species"))

write.csv(outdata,file=outfile,quote=FALSE, row.names=FALSE)

