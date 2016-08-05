#ifndef MESH_H
#define MESH_H

void cube_cartesian(vector <voxel>& voxels,simulation *sys,double D,int M){


    double h = pow(sys->volume,1.0/3.0)/M;
    sys->volume = 0.0;
    voxel temp;
    temp.vol = pow(h,3);
    temp.sd = 0;
    temp.totalD = 6.0*D/pow(h,2.0);
    temp.neighbors.resize(6);
    temp.num_conn = 6;
    temp.neighbors[0].D = D/pow(h,2.0);
    temp.neighbors[1].D = D/pow(h,2.0);
    temp.neighbors[2].D = D/pow(h,2.0);
    temp.neighbors[3].D = D/pow(h,2.0);
    temp.neighbors[4].D = D/pow(h,2.0);
    temp.neighbors[5].D = D/pow(h,2.0);
    
    
    
    for(int k=0;k<M;k++){
        for(int j=0;j<M;j++){
            for(int i=0;i<M;i++){

                temp.id = i+M*j+M*M*k;
                temp.node[0] = (i+0.5)*h;
                temp.node[1] = (j+0.5)*h;
                temp.node[2] = (k+0.5)*h;

                sys->volume += temp.vol;
                
                temp.neighbors[0].vox = (i+1)+M*j+M*M*k;
                temp.neighbors[1].vox = (i-1)+M*j+M*M*k;
                temp.neighbors[2].vox = i+M*(j+1)+M*M*k;
                temp.neighbors[3].vox = i+M*(j-1)+M*M*k;
                temp.neighbors[4].vox = i+M*j+M*M*(k+1);
                temp.neighbors[5].vox = i+M*j+M*M*(k-1);

                if(i==0){
                    temp.neighbors[1].vox += 1;
                }
                else if(i==M-1){
                    temp.neighbors[0].vox -= 1;
                }

                if(j==0){
                    temp.neighbors[3].vox += M;
                }
                else if(j==M-1){
                    temp.neighbors[2].vox -= M;
                }

                if(k==0){
                    temp.neighbors[5].vox += M*M;
                }
                else if(k==M-1){
                    temp.neighbors[4].vox -= M*M;
                }

                voxels.push_back(temp);

            }
        }
    }
    printf("system volume = %g\n",sys->volume);
}

void read_mesh(FILE *fin,vector <voxel>& voxels,simulation *sys){
    char line[MAX_CHARACTERS];
    memset(line,'\0',sizeof(line));

    char *ptr;
    neighbor tempn;
    double total_vol = 0.0;
    while(!feof(fin)){
        voxel temp;
        if(fgets(line,MAX_CHARACTERS,fin)!=NULL){
//                        printf("%s\n",line);
//            printf("1\n");
            //        ti = strtok(line," ");
            temp.id = atoi(strtok(line," "))-1;
//            printf("2\n");
            //            printf("%d\n",temp.id);
            temp.node[0] = strtod(strtok(NULL," "),&ptr);
            temp.node[1] = strtod(strtok(NULL," "),&ptr);
            temp.node[2] = strtod(strtok(NULL," "),&ptr);
//            printf("3\n");
            temp.num_conn = atoi(strtok(NULL," "));
            temp.totalD = 0.0;
//            printf("4\n");
            for(int i=0;i<temp.num_conn;i++){
                tempn.vox = atoi(strtok(NULL," "))-1;
                tempn.D = strtod(strtok(NULL," "),&ptr);
                temp.neighbors.push_back(tempn);
                //                printf("tempn.D=%g\n",tempn.D);
                temp.totalD += tempn.D;
            }
//            printf("5\n");
            //            printf("*********\n");
            //            printf("Voxel =%d\n",temp.id);
            //            printf("totalD=%g\n",temp.totalD);
            //            printf("volume=%g\n",temp.vol);
            temp.vol = strtod(strtok(NULL," "),&ptr);
            total_vol += temp.vol;
            temp.sd = atoi(strtok(NULL," "));
//            printf("2\n");
            voxels.push_back(temp);
//            printf("3\n");
        }
    }
    sys->volume = total_vol;
}

#endif
