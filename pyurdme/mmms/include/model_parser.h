#ifndef MODEL_PARSER_H
#define MODEL_PARSER_H

#include "stdlib.h"

void print_model(vector <species>& specs,vector <association>& associations,vector <dissociation>& dissociations,vector <birth>& births,simulation *sim){


    printf("********************\n");
    printf("Number of species: %d\n",(int)(specs.size()));
    for(int i=0;i<(int)(specs.size());i++){
        printf("--------------------\n");
        printf("Species %s:\n",specs[i].name);
        printf("D=%.5g\n",specs[i].D);
        printf("sigma=%.5g\n",specs[i].sigma);
        printf("Dimension=%d\n",specs[i].dim);
        printf("Initial value: %d\n",specs[i].initial_value);
    }
    printf("********************\n");
    printf("Number of reactions: %d\n",(int)(associations.size())+(int)(dissociations.size()+(int)(births.size())));
    printf("--------------------\n");
    printf("Bimolecular:\n");
    printf("--------------------\n");
    for(int i=0;i<(int)(associations.size());i++){
        printf("%s+%s-->",specs[associations[i].reactant1].name,specs[associations[i].reactant2].name);
        if((int)(associations[i].products.size())==1){
            printf("%s, k=%.5g\n",specs[associations[i].products[0]].name,associations[i].k);
        }
        else if((int)(associations[i].products.size())==2){
            printf("S%d+S%d, k=%.5g\n",associations[i].products[0],associations[i].products[1],associations[i].k);
        }
    }
    printf("--------------------\n");
    printf("Unimolecular:\n");
    printf("--------------------\n");
    for(int i=0;i<(int)(dissociations.size());i++){
        printf("%s-->",specs[dissociations[i].reactant].name);
        if((int)(dissociations[i].products.size())==1){
            printf("%s, k=%.5g\n",specs[dissociations[i].products[0]].name,dissociations[i].k);
        }
        else if((int)(dissociations[i].products.size())==2){
            printf("%s+%s, k=%.5g\n",specs[dissociations[i].products[0]].name,specs[dissociations[i].products[1]].name,dissociations[i].k);
        }
        else if((int)(dissociations[i].products.size())==0){
            printf("0, k=%.5g\n",dissociations[i].k);
        }
    }
    printf("--------------------\n");
    printf("Birth processes:\n");
    printf("--------------------\n");
    for(int i=0;i<(int)(births.size());i++){
        printf("0-->");
        printf("%s, k=%.5g\n",specs[births[i].product].name,births[i].k);

    }
    printf("********************\n");
    printf("Volume: %.5g\n",sim->volume);
    printf("Boundary: [%g %g %g %g",sim->boundary[0],sim->boundary[1],sim->boundary[2],sim->boundary[3]);
    if(sim->dimension==3){
        printf(" %g %g]\n",sim->boundary[4],sim->boundary[5]);
    }
    else{
        printf("]\n");
    }
    printf("T: %.5g\n",sim->T);
    printf("Number of intervals: %d\n",sim->num_intervals);


}


void print_species(vector <species>& specs){
    for(vector<species>::iterator it = specs.begin(); it != specs.end(); ++it) {
        printf("Species name: %s\n",it->name);
        printf("---Diffusion constant: %g\n",it->D);
        printf("---Reaction radius: %g\n",it->sigma);
        printf("---Initial value: %d\n",it->initial_value);
    }
}

void print_parameters(vector <parameter>& params){
    for(vector<parameter>::iterator it = params.begin(); it != params.end(); ++it) {
        printf("Parameter name: %s\n",it->name);
        printf("---Value: %g\n",it->value);
    }
}

void print_reactions(vector <association>& assocs,vector <dissociation>& dissocs,vector <birth>& births){

    for(int i=0;i<(int)(assocs.size());i++){
        printf("%d+%d->",assocs[i].reactant1,assocs[i].reactant2);
        for(int j=0;j<(int)(assocs[i].products.size());j++){
            printf("%d",assocs[i].products[j]);
            if(j<(int)(assocs[i].products.size())-1){
                printf("+");
            }
        }
        printf(", k=%g\n",assocs[i].k);
    }
    for(int i=0;i<(int)(dissocs.size());i++){
        printf("%d->",dissocs[i].reactant);
        for(int j=0;j<(int)(dissocs[i].products.size());j++){
            printf("%d",dissocs[i].products[j]);
            if(j<(int)(dissocs[i].products.size())-1){
                printf("+");
            }
        }
        printf(", k=%g\n",dissocs[i].k);
    }
    for(int i=0;i<(int)(births.size());i++){
        printf("->%d, k=%g\n",births[i].product,births[i].k);
    }

}

/* Attempts to parse 'input', which can be either a numeric value or a predefined parameter.
 * On failure: status =-1.
 * On success: status = 1.
 */
double parse_parameter(char *input,vector <parameter>& params,int *status){
    char *ptr;
//    printf("INPUT %s\n",input);
    double test = strtod(input,&ptr);
    if(test!=0){
        *status = 1;
        return test;
    }
    if(test==0 && (strcmp(input,"0")==0 || strcmp(input,"0.0")==0)){
        *status = 1;
        return test;
    }
    for(vector<parameter>::iterator it = params.begin(); it != params.end(); ++it){
        if(strcmp(input,it->name)==0){
            *status = 1;
            return it->value;
        }
    }
    *status = -1;
    return 0.0;
}

/* Attempts to find species 'spec'. Returns index upon success, -1 on failure. */
int parse_species(char *spec,vector <species>& specs){
    int N = (int)(specs.size());
    for(int i=0;i<N;i++){
        if(strcmp(spec,specs[i].name)==0){
            return i;
        }
    }
    printf("Failed to parse reactions. Could not find species matching: %s\n",spec);
    exit(EXIT_FAILURE);
}


void add_reaction(char *sp,vector <association>& assocs,vector <dissociation>& dissocs,vector <birth>& births,vector <species>& specs,vector <parameter>& params){

    printf("reading reaction (1): %s\n",sp);
    char *r;
    r = strtok(sp,">");
    printf("reading reaction (1): %s\n",r);
    if(r==NULL){
        printf("Failed to parse reaction: %s\n",sp);
        exit(EXIT_FAILURE);
    }

    char *p;
    p = strtok(NULL,">");
    printf("reading reaction (2): %s\n",p);
    if(p==NULL){
        p = r;
        r = NULL;
    }


    vector <char*> reactants;
    vector <char*> products;

    double rrate;

    /* Parse reactants. */

    r = strtok(r," ");
    while(r!=NULL){
        reactants.push_back(r);
        r = strtok(NULL," ");
    }
    int num_reactants = (int)(reactants.size());

    /* Parse products. */
    p = strtok(p," ");
    while(p!=NULL){
        products.push_back(p);
        p = strtok(NULL," ");
    }
    int num_products = (int)(products.size());
    printf("reading reaction: %s, %s\n",r,p);
    /* Parse reaction rate. */
    int status;
    rrate = parse_parameter(products[num_products-1],params,&status);
    if(status==-1 || rrate<0){
        printf("Failed to read reaction rate: %s\n",products[num_products-1]);
        exit(EXIT_FAILURE);
    }

    if(num_reactants==2){
        association temp_reac;
        temp_reac.reactant1 = parse_species(reactants[0],specs);
        temp_reac.reactant2 = parse_species(reactants[1],specs);
        for(int i=0;i<num_products-1;i++){
            temp_reac.products.push_back(parse_species(products[i],specs));
        }
        temp_reac.k = rrate;
        assocs.push_back(temp_reac);
    }
    else if(num_reactants==1){
        dissociation temp_reac;
        temp_reac.reactant = parse_species(reactants[0],specs);
        for(int i=0;i<num_products-1;i++){
            temp_reac.products.push_back(parse_species(products[i],specs));
        }
        temp_reac.k = rrate;
        dissocs.push_back(temp_reac);
    }
    else if(num_reactants==0){
        birth temp_reac;
        if(num_products-1!=1){
            printf("A birth process can only have one product.\n");
            exit(EXIT_FAILURE);
        }
        temp_reac.product = parse_species(products[0],specs);
        temp_reac.k = rrate;
        births.push_back(temp_reac);
    }

}

void add_parameter(char *sp,vector <parameter>& params){
    char *ptr;
    char *str;
    str = strtok(sp," ");
    if(str==NULL){
        printf("Error when reading data for parameter: %s (1)\n",sp);
        exit(EXIT_FAILURE);
    }

    parameter temp;
    memset(temp.name,'\0',sizeof(temp.name));
    strcpy(temp.name,str);
    str = strtok(NULL," ");
    if(str==NULL){
        printf("Error when reading parameter: %s (2)\n",temp.name);
        exit(EXIT_FAILURE);
    }
    temp.value = strtod(str,&ptr);
    params.push_back(temp);
}

void add_species(char *sp,vector <species>& specs,vector <parameter>& params){
    char *ptr;
    char *str;
    str = strtok(sp," ");
    if(str==NULL){
        printf("Error when reading data for species: %s (1)\n",sp);
        exit(EXIT_FAILURE);
    }

    species temp;
    memset(temp.name,'\0',sizeof(temp.name));
    strcpy(temp.name,str);

    /* TODO: Code is being repeated here. */

    str = strtok(NULL," ");
    if(str==NULL){
        printf("Error when reading data for species: %s (2)\n",temp.name);
        exit(EXIT_FAILURE);
    }
    int status;
    temp.D = parse_parameter(str,params,&status);
    if(status==-1 || temp.D<0){
        printf("Invalid diffusion constant for species %s: %s\n",temp.name,str);
        exit(EXIT_FAILURE);
    }

    str = strtok(NULL," ");
    if(str==NULL){
        printf("Error when reading data for species: %s (3)\n",sp);
        exit(EXIT_FAILURE);
    }
    temp.sigma = parse_parameter(str,params,&status);
    if(status==-1 || temp.sigma<=0){
        printf("Invalid reaction radius for species %s: %s\n",temp.name,str);
        exit(EXIT_FAILURE);
    }

    str = strtok(NULL," ");
    if(str==NULL){
        printf("Error when reading data for species: %s (4)\n",sp);
        exit(EXIT_FAILURE);
    }
    temp.initial_value = strtol(str,&ptr,10);
    if(temp.initial_value<0){
        printf("Invalid initial value for species %s: %s\n",temp.name,str);
        exit(EXIT_FAILURE);
    }
    str = strtok(NULL," ");
    if(str==NULL){
        printf("Error when reading data for species: %s (4)\n",sp);
        exit(EXIT_FAILURE);
    }
    temp.meso_micro = strtol(str,&ptr,10);
    if(temp.meso_micro!=0 && temp.meso_micro!=1){
        printf("Invalid modeling level for species %s: %s\n",temp.name,str);
        exit(EXIT_FAILURE);
    }
    
    specs.push_back(temp);
    printf("Added species: %s\n",temp.name);
}

char *fetch_next(char *input,char *err_msg,bool first){
    char *temp;
    if(first){
        temp = strtok(input," ");
    }
    else{
        temp = strtok(NULL," ");
    }
//    printf("Fetched: %s\n",temp);
//    str = &temp;
    if(temp==NULL){
        printf("%s: %s\n",err_msg,input);
        exit(EXIT_FAILURE);
    }
    return temp;
}

void add_boundary(char *input,simulation *sim,vector <parameter>& params){
    char *str;



    int status;
    char err_msg[MAX_CHARACTERS];

    for(int i=0;i<2*sim->dimension;i++){
        sprintf(err_msg,"Error reading boundary (%d)",i);
        if(i==0)
            str = fetch_next(input,err_msg,true);
        else
            str = fetch_next(input,err_msg,false);
//        printf("To parse parameter: %s\n",str);
        sim->boundary[i] = parse_parameter(str,params,&status);
    }

}

/* Function that parses model and sets up basic structures. See documentation for
 * specification of model format.
 * Input:
 * Output:
 */
void parse_model(char *filename,simulation *sim,vector <species>& specs,vector <association>& assocs,vector <dissociation>& dissocs,vector <birth>& births,vector <parameter>& parameters){

    const int NUM_KEYWORDS = 11;
    const char *keywords[] = {"NAME","SPECIES","PARAMETER","REACTION","DIMENSION","BOUNDARY","OUTPUT","NTRAJ","NCORES","T","NINT"};

    FILE *model_in;
    model_in = fopen(filename,"r");

    bool dimension_set = false;

    if(model_in==NULL){
        printf("Failed to open model file.\n");
        exit(EXIT_FAILURE);
    }
    else{
        printf("Parsing model...\n");

        char line[MAX_CHARACTERS];
        line[MAX_CHARACTERS-2] = '\0';

        while(!feof(model_in)){
            if(fgets(line,MAX_CHARACTERS,model_in)==NULL){
                 if(feof(model_in))
                     goto end;
                 printf("Error reading model. Exiting.\n");
                 printf("%s\n",line);
                 printf("NINT=%d\n",sim->num_intervals);
                 exit(EXIT_FAILURE);
             }

//            printf("%s\n",line);
//            fgets(line,MAX_CHARACTERS,model_in);

            /* Ignore blank lines and comments. */
            while(line[0]=='\0' || line[0]=='\n' || line[0]=='#'){
                if(feof(model_in))
                    goto end;


                if(fgets(line,MAX_CHARACTERS,model_in)==NULL){
                    if(feof(model_in))
                        goto end;
                    printf("Error reading model. Exiting. 2.\n");
                    printf("%s\n",line);
                    exit(EXIT_FAILURE);
                }

            }

            /* Make sure that the line is less than MAX_CHARACTERS characters long. */
            if(line[MAX_CHARACTERS-2]!='\0'){
                printf("Error reading model file. Line exceeds maximum length of %d characters.\nNear: %s\n",MAX_CHARACTERS-1,line);
                exit(EXIT_FAILURE);
            }


            char *keyword;
            keyword = strtok(line," ");
            char *input;
            input = strtok(NULL,"\n");

            int keyword_loc = -1;
            char *ptr;

            /* Search for the command. */
            for(int i=0;i<NUM_KEYWORDS;i++){
                if(strcmp(keyword,keywords[i])==0){
//                    printf("Reading: %s...\n",keyword);
//                    printf("Input is: %s\n",input);
                    keyword_loc = i;
                    break;
                }
            }
            if(keyword_loc==-1){
                printf("Unrecognized command near\n%s\n",keyword);
                exit(EXIT_FAILURE);
            }
            else if(strcmp(keywords[keyword_loc],"NAME")==0){
                memset(sim->name,'\0',sizeof(sim->name));
                strcpy(sim->name,input);
                printf("Name of model: %s\n",sim->name);
            }
            else if(strcmp(keywords[keyword_loc],"NTRAJ")==0){
//                printf("%s\n",input);
                sim->ntraj = strtol(input,&ptr,10);
                if(sim->ntraj>0){
                    printf("Number of trajectories: %d\n",sim->ntraj);
                }
                else{
                    printf("Incorrect number of trajectories specified:\n%s\n",input);
                    exit(EXIT_FAILURE);
                }
            }
            else if(strcmp(keywords[keyword_loc],"NCORES")==0){
                sim->ncores = strtol(input,&ptr,10);
                if(sim->ncores>0){
                    printf("Number of cores: %d\n",sim->ncores);
                }
                else{
                    printf("Incorrect number of cores specified:\n%s\n",input);
                    exit(EXIT_FAILURE);
                }
            }
            else if(strcmp(keywords[keyword_loc],"DIMENSION")==0){
                sim->dimension = strtol(input,&ptr,10);
                if(sim->dimension==2 || sim->dimension==3){
                    printf("Dimension: %d\n",sim->dimension);
                }
                else{
                    printf("Incorrect dimension specified:\n%s\n",input);
                    exit(EXIT_FAILURE);
                }
                dimension_set = true;
            }
            else if(strcmp(keywords[keyword_loc],"SPECIES")==0){
                add_species(input,specs,parameters);
            }
            else if(strcmp(keywords[keyword_loc],"PARAMETER")==0){
                add_parameter(input,parameters);
            }
            else if(strcmp(keywords[keyword_loc],"REACTION")==0){
                add_reaction(input,assocs,dissocs,births,specs,parameters);
            }
            else if(strcmp(keywords[keyword_loc],"T")==0){
                sim->T = strtod(input,&ptr);
                if(sim->T<=0){
                    printf("ERROR: Final time not valid: %s\n",input);
                    exit(EXIT_FAILURE);
                }
            }
            else if(strcmp(keywords[keyword_loc],"NINT")==0){
                sim->num_intervals = strtol(input,&ptr,10);
                if(sim->num_intervals<=0){
                    printf("ERROR: Number of intervals not valid: %s\n",input);
                    exit(EXIT_FAILURE);
                }
            }
            else if(strcmp(keywords[keyword_loc],"BOUNDARY")==0){
                if(!dimension_set){
                    printf("ERROR: Trying to set boundary, but dimension undefined.\n");
                    exit(EXIT_FAILURE);
                }
                add_boundary(input,sim,parameters);
                sim->volume = 1.0;
                for(int i=0;i<sim->dimension;i++){
                    sim->volume *= sim->boundary[2*i+1]-sim->boundary[2*i];
                }
            }
            else if(strcmp(keywords[keyword_loc],"OUTPUT")==0){
                printf("reading output settings...\n");
            }
        }
        end:
        print_species(specs);
        print_parameters(parameters);
        print_reactions(assocs,dissocs,births);
        printf("Finished reading model file.\n");

    }
}


#endif
