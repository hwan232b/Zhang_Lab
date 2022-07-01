#create 2 matrix
m1 <- matrix(1:9, nrow=3, ncol=3, dimnames = list(c("r1","r2","r3"), c("c2","c1","c3")))
print(m1)

m2 <- matrix(1:15, nrow=5, ncol=3, dimnames = list(c("r1","r2","r3","r4","r5"), c("cell","c2","c3")))
print(m2)

m3 <- matrix(nrow=3, ncol=3, dimnames = list(c("r1",NA,NA), c("cell",NA,"c3")))
print(m3)

## set the value for the first column and every row to Cn (c1, c2, ...)
k <- 1
for(i in 4:1){
  for(j in 1:1){
    t <- paste0("c",i)
    print(t)
    m2[k,j] <- t
  }
  k <- k+1
}

#conver matrix to data frame
m1_df <- as.data.frame(m1)
m2_df <- as.data.frame(m2)
m3_df <- as.data.frame(m3)

intersect2 <- intersect(colnames(m1_df),m2_df$cell)
m3_df <- m2_df[m2_df$cell %in% intersect2,]

#compare, expect to be FALSE
all(colnames(m1_df) == m3_df$cell)

m1_df <- m1_df[, order(colnames(m1_df))]

#expect to be TRUE
all(colnames(m1_df) == m3_dfaz$cell)

m3_df[is.na(x = m3_df)] <- 0
df <-data.frame(m1_df,m3_df)
NA %in% m3