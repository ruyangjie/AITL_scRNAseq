

#' @title 生成差异较大的颜色系列
#' @description：
#' 修饰过的rainbow函数，用于生成差异较大、数量较多的颜色系列
#'
#' @param    n 指定色阶的个数
#' @return   一个包含颜色信息的向量（vector）
#' @export


dif_seq_rainbow<-function(n){

  library("colorspace")

  #n=9
  div=11
  depth= n %/% div+2
  if (n %% div!=0) depth=depth+1

  base_color<-c("red","#FF6B00","gold","greenyellow","green","cyan","#007FFF","blue","#4A00E9","darkviolet","#C90069")
  color_matrix<-matrix(0,ncol=depth,nrow=0)

  for(col_i in base_color){

    seq_col<-colorRampPalette(c("white",col_i,"black"))(depth)

    color_matrix<-rbind(color_matrix,seq_col)
  }

  color_matrix<-color_matrix[,c(-1,-1*depth)]

  color_matrix_order<-data.frame(n=c(1:(depth-2)),t(color_matrix))
  color_matrix_order$n<-abs(color_matrix_order$n-mean(color_matrix_order$n))
  color_matrix_order2<-color_matrix_order[order(color_matrix_order[,"n"]),]
  color_matrix_order3<-t(color_matrix_order2[,-1])

  set.seed(456)
  rdm<-function(var){sample(var,length(var))}
  color_matrix_order4<-apply(color_matrix_order3,2,rdm)

  color_vector<-as.vector(color_matrix_order4)[c(1:n)]
  #if(n>30) cat("NOTE:You are using dif_seq_hcl to produce more than 30 colors, some of them may be difficult to distinguish\n")
  return(color_vector)
}


