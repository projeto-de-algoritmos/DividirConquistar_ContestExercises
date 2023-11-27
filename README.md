# Exercícios Online

**Conteúdo da Disciplina**: Dividir e Conquistar <br>

## Alunos
|Matrícula | Aluno |
| -- | -- |
| 21/1030729  |  Eric Rabelo Borges |
| 19/0058455  |  Mariana Oliveira Pires do Rio |

## Sobre 

Este projeto abrange as soluções para exercícios de Juiz Online que demandam a aplicação de técnicas de dividir e conquistar.
 
### Exercício 1:  [Golf Bot](http://uva.onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&category=24&page=show_problem&problem=4744)
O exercício aborda a perspectiva de alguém que, inicialmente, expressa aversão ao golfe e decide criar um robô de golfe infalível. Contudo, o robô enfrenta desafios durante o jogo devido à limitação do botão seletor de distância embutido, que não contempla todas as opções necessárias. Agora, no campo de golfe, a pessoa busca determinar quantos buracos o robô pode completar em duas tacadas ou menos, enfrentando o desafio de otimizar o desempenho do robô diante das restrições, especialmente a falta de opções adequadas de distância.
O exercicio foi resolvido usando FFT.

Confira o enunciado do exercício [aqui](http://uva.onlinejudge.org/index.php?option=com_onlinejudge&Itemid=8&category=24&page=show_problem&problem=4744)

- [Ir para solução](Solutions/golf_bot.cpp)

### Exercício 2: [Maximum Self-Matching](https://www.spoj.com/problems/MAXMATCH/)

O exercício consiste em uma string s composta por letras 'a', 'b' e 'c'. A função de correspondência ms(i) é definida como o número de caracteres coincidentes entre a string original e sua cópia deslocada por i posições. A ideia consiste em calcular o valor máximo de ms(i) para todas as posições i no intervalo de 1 até o comprimento da string e apresentar todas as posições ótimas de i em ordem crescente. 
A solução demanda identificar as posições que maximizam o número de caracteres correspondentes na string e ordená-las para fornecer uma resposta completa ao exercício proposto, o exercicio foi resolvido usando FFT.

Confira o enunciado do exercício [aqui](https://www.spoj.com/problems/MAXMATCH/)

- [Ir para solução](Solutions/maxmatch.cpp)
### Exercício 3: [Dr. Evil Underscores](https://codeforces.com/problemset/problem/1285/D)

O problema consiste em encontrar o valor mínimo possível da operação XOR bit a bit entre um inteiro X e um conjunto de inteiros dados. Mais especificamente, dado um conjunto de n inteiros, o objetivo é escolher um inteiro X de forma a minimizar o valor máximo da operação XOR entre X e qualquer elemento do conjunto. O código C++ fornecido utiliza uma abordagem de dividir e conquistar para explorar recursivamente possíveis valores de X, dividindo o conjunto em subconjuntos com base na posição do bit mais significativo.

Confira o enunciado do exercício [aqui](https://codeforces.com/problemset/problem/1285/D)

- [Ir para solução](Solutions/Dr.EvilUnderscores.cpp)

### Exercício 4: [a-Good String](https://codeforces.com/contest/1385/problem/D)

A questão consiste em determinar o número mínimo de movimentos necessários para transformar uma string dada em uma "a-good string", onde uma "a-good string" atende a condições específicas relacionadas à presença do caractere 'a' e à sua relação com os demais caracteres. A solução proposta utiliza uma abordagem de dividir e conquistar, onde a função menorCusto é recursivamente chamada para calcular o custo mínimo de movimentos em diferentes intervalos da string, considerando as condições de formação da "a-good string". 

Confira o enunciado do exercício [aqui](https://codeforces.com/contest/1385/problem/D)

- [Ir para solução](Solutions/a-GoodString.cpp)


## Screenshots

### Golf Bot
![](Assets/Golf_Bot.png)

### Maximum Self-Matching
![](Assets/max.png)

### Dr. Evil Underscores
![](Assets/Dr.EvilUnderscores.png)

### a-Good String
![](Assets/aGoodString.cpp.png)

## Instalação 
Não é necessario instalação. As soluções podem ser submetidas nos sites: [Codeforces](https://codeforces.com/), [UVa](https://onlinejudge.org/index.php) e [Spoj](https://www.spoj.com/).

**Linguagem**: C++<br>

## Apresentação
Assista o vídeo explicativo do projeto clicando [aqui](/Apresentação.mp4). <br>
*Observação: A apresentação está em formato .mp4, portanto, faça o download para assisti-la.*