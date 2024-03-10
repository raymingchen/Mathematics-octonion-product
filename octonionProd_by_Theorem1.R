#### x_8 octonion product (in Theorem 1) on causlity######
OctProd=function(vecW,vecV,n) ##|vecW|=|vecV|=8;n=n'th coefficient of the product
{
vecW_n_plus=c(vecW[0+1],vecW[(n+4)%%8+1],vecW[(n+1)%%8+1],vecW[(n+2)%%8+1]);
vecW_n_minus=c(-vecW[n+1],vecW[(n+5)%%8+1],vecW[(n+3)%%8+1],vecW[(n+6)%%8+1]);
vecV_n_plus=c(vecV[n+1],vecV[(n+5)%%8+1],vecV[(n+3)%%8+1],vecV[(n+6)%%8+1]);
vecV_n_minus=c(vecV[0+1],vecV[(n+4)%%8+1],vecV[(n+1)%%8+1],vecV[(n+2)%%8+1]);
res_3terms=vecW_n_plus*vecV_n_plus-vecW_n_minus*vecV_n_minus;
return(sum(res_3terms));
}

################# data analysis #################
vecW=c(180,23,28,36,53,67,88,92); 
vecW2=c(168,9,20,36,52,79,100,117);
vecV=c(200,12,45,47,62,88,110,157); 
vecV2=c(189,21,32,40,55,65,81,132);
cof_0=vecW[1]*vecV[1]-vecW[2:8]%*%vecV[2:8];
cof2_0=vecW2[1]*vecV2[1]-vecW2[2:8]%*%vecV2[2:8];

newVEC=rep(0,8);
newVEC2=rep(0,8);
newVEC[1]=cof_0;
newVEC2[1]=cof2_0; 
for(n in 1:7)
{
newVEC[n+1]=OctProd(vecW,vecV,n);
newVEC2[n+1]=OctProd(vecW2,vecV2,n);
}
den=sum(newVEC[2:8]-newVEC[1]); 
den2=sum(newVEC2[2:8]-newVEC2[1]);
weight=(newVEC[2:8]-newVEC[1])/den;
weight2=(newVEC2[2:8]-newVEC2[1])/den2;
causality=weight%*%newVEC[2:8];
causality2=weight2%*%newVEC2[2:8];
res1=causality/abs(newVEC[1]);
res2=causality2/abs(newVEC2[1]);

res1;
res2;


