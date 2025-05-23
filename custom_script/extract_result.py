import os 
from sys import argv

#在CTDG_out 目录下运行,会将spe_list中所有物种的结果提取出来存放到一个文件中
spe_list = argv[1]
with open(spe_list,'r') as f:
    spe_name = [line.strip() for line in f]
work_path = os.getcwd() 
print("species:\n"+str(spe_name))
 
#读取genes_clean文件
with open("genes_result.txt","w") as out:  
    for spe in spe_name:   
        # print(spe)    
        if os.path.isdir(spe):
            spe_path=os.path.join(work_path,spe)   
            gene_familys = os.listdir(spe_path)
            for family in gene_familys:
                family_path=os.path.join(spe_path,family)
                report_path=os.path.join(family_path,"report")
                numbers_file_name=family+"_genes_clean.csv" 
                numbers_file_path=os.path.join(report_path,numbers_file_name)
                if os.path.isfile(numbers_file_path):
                    #读取结果文件        
                    with open(numbers_file_path,"r") as infile:
                        infile.readline()
                        for line in infile:
                            acc = line.split(",")[0]
                            species = line.split(",")[1]
                            chromosome = line.split(",")[2]
                            cluster = line.split(",")[3]
                            order = line.split(",")[-1]
                            print("%s\t%s\t%s\t%s\t%s"%(acc,species,chromosome,cluster,order),end = "",file = out)
        else:
            continue
       

## 读取numbers文件
with open("numbers_result.txt",'w') as out:  
    for spe in spe_name:   
        # print(spe)    
        if os.path.isdir(spe):
            spe_path=os.path.join(work_path,spe)   
            gene_familys = os.listdir(spe_path)
            for family in gene_familys:
                family_path=os.path.join(spe_path,family)
                report_path=os.path.join(family_path,"report")
                numbers_file_name=family+"_numbers_clean.csv" 
                numbers_file_path=os.path.join(report_path,numbers_file_name)
                if os.path.isfile(numbers_file_path):
                    #读取结果文件        
                    with open(numbers_file_path,"r") as infile:
                        infile.readline()
                        for line in infile:
                            species = line.split(",")[0]
                            chromosome=line.split(",")[1]
                            cluster=line.split(",")[2]
                            gene_duplicates=line.split(",")[3]
                            start=line.split(",")[4]
                            end=line.split(",")[5]
                            length=line.split(",")[6]
                            p_95=line.split(",")[7]
                            print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(species,family,chromosome,cluster,gene_duplicates,start,end,length,p_95),end="",file=out)
                            

        else:
            continue                    
                    
                    





