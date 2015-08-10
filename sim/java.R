library(rJava)
.jinit()

java.version <- function(){
  print(.jcall("java/lang/System", "S", "getProperty", "java.runtime.version"))
}