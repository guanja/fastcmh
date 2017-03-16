#script to create sample data
createSampleData <- function(){
    require(fastcmh)

    makefastcmhdata(folder="./", L=20, n=50, K=2, tau1=5, taulength1=4, tau2=12, taulength2=4, seednum=3)

    cat("done\n")
}
