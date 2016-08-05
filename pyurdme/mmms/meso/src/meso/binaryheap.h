#ifndef BINARYHEAP
#define BINARYHEAP

#include <vector>

using namespace std;


//typedef struct tent_event{
//    int type; //-1 diffusion, 0 dissociation, 1 association
//    int index;
//    double t;
//    vector <int> reactants;
//    bool active;
//}tent_event;


//Note that second child is simply i+1 if first child is i.
int first_child(int index){
    return 2*index+1;
}

int second_child(int index){
    return 2*index+2;
}

int parent(int index){
    return floor(((double)index-1.0)/2.0);
}

void heap_swap_elems(vector <tent_event>& elems,int index1,int index2){
    iter_swap(elems.begin()+index1,elems.begin()+index2);
}

void heap_move_up(vector <tent_event>& elems,int index){
    heap_swap_elems(elems,index,parent(index));
}

void heap_sift_down(vector <tent_event>& elems,int index,bool (*compare)(tent_event*,tent_event*)){
    int N = (int)(elems.size());
    int c1,c2,minindex;
    c1 = first_child(index);
    c2 = second_child(index);
    if(c2>N-1){
        if(c1>N-1){
            return;
        }
        else{
            minindex = c1;
        }
    }
    else{
        if(compare(&elems[c1],&elems[c2])){
            minindex = c1;
        }
        else{
            minindex = c2;
        }
    }
    if(compare(&elems[minindex],&elems[index])){
        heap_swap_elems(elems,index,minindex);
        heap_sift_down(elems,minindex,compare);
    }
    
}

//void heap_del_head(vector <tent_event>& elems,bool (*compare)(tent_event*,tent_event*)){
//    int N = (int)(elems.size());
////    heap_swap_elems(elems,0,N-1);
//    elems[0] = elems[N-1];
//    elems.erase(elems.end());
//    heap_sift_down(elems,0,compare);
//    
//}


//void initialize_heap(vector <tent_event>& elems){
//    
//}

void heap_insert(vector <tent_event>& elems,tent_event new_elem,bool (*compare)(tent_event*,tent_event*)){
    elems.push_back(new_elem);
    int N = (int)(elems.size());
    int index = N-1;
    while(index>0 &&compare(&new_elem,&elems[parent(index)])){
        heap_move_up(elems,index);
        index = parent(index);
    }
}

void heap_delete(vector <tent_event>& elems,int index,bool (*compare)(tent_event*,tent_event*)){
    int N = (int)(elems.size());
    elems[index] = elems[N-1];
    elems.erase(elems.end());
    heap_sift_down(elems,index,compare);
}

void print_heap(vector <tent_event>& elems){
    int N = (int)(elems.size());
    for(int i=0;i<N;i++){
        printf("%g\n",elems[i].t);
    }
}

void unroll_heap(vector <tent_event>& elems,bool (*compare)(tent_event*,tent_event*)){
    while((int)(elems.size())>0){
        printf("%g\n",elems[0].t);
        heap_delete(elems,0,compare);
    }
    
}


#endif
