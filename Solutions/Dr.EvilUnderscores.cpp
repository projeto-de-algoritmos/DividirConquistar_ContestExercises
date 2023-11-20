#include <bits/stdc++.h>
using namespace std;

long long busca(set<long long> &a, int posicao){
    //ja atingimos o final do numero
    if(posicao == 0)
        return 0; 

    //preenchemos os sets 0 e 1 
    set<long long> set0, set1;  
    long long mask = 1LL << (posicao - 1); 
    for(long long i : a) {
        if((mask & i) == 0)
            set0.insert(i); 
        else 
            set1.insert(i);
    }

    //se o set 0 estiver vazio, chamamos a funcao recursiva para o set 1
    if(set0.size() == 0)
        return busca(set1, posicao - 1);

    //se o set 1 estiver vazio, chamamos a funcao recursiva para o set 0
    if(set1.size() == 0)
        return busca(set0, posicao - 1);
    
    //caso nenhum dos sets esteja vazio, retornamos o minimo entre o set 0 e 1
    return (1LL << (posicao - 1)) + min(busca(set0, posicao - 1),
                                        busca(set1, posicao - 1));
}

int main(){
    int n; 
    cin >> n; 

    set<long long> a; 
    for(int i = 0; i < n; i++){
        long long t; 
        cin >> t; 
        a.insert(t); 
    }

    cout << busca(a, 30); 
    return 0; 
}