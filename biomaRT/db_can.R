library

db<-read.csv("/Users/nicholaswaters/Documents/GU_R/biomaRT/CAdb.csv")
for (i in row.names(db)){
  if (length(grep(pattern = "orf19.*", x = db[i,], ignore.case = TRUE, perl = TRUE, value = TRUE))==1){
    db$orf19a<-grep(pattern = "orf19.*", x = row, ignore.case = TRUE, perl = TRUE, value = TRUE)
  }
  else{
    print("error: you've been duped by duplicates!")}
}
