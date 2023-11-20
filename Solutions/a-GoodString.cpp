#include <bits/stdc++.h>
using namespace std;
char caracteres[150000]; 

//checa quantos caracteres precisam ser trocados
int custo(int inicio, int fim, char c){
    int custo = 0;
    for(int i = inicio; i <= fim; i++){
        if(caracteres[i] != c)
            custo++;
    }
    return custo;
}

int menorCusto(int inicio, int fim, char c){
    //se o intervalo tiver tamanho 1, checa se o caractere é o desejado
    if(inicio == fim){
        if(caracteres[inicio] == c)
            return 0;
        else
            return 1;
    }

    int meio = (inicio + fim) / 2;
    
    return min(custo(inicio, meio, c) + menorCusto(meio + 1, fim, c + 1),
                custo(meio + 1, fim, c) + menorCusto(inicio, meio, c + 1));
}

int main(){
    int t; 
    cin >> t;

    while(t--){
        int n; 
        cin >> n;

        for(int i = 1; i <= n; i++)
            cin >> caracteres[i];

        cout << menorCusto(1, n, 'a') << endl;
    }
    
    return 0; 
}