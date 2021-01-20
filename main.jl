using LinearAlgebra
# 10x + 3y - 2z = 57
# 2x + 8y -1z = 20
# x + y + 5z = -4
# chute inicial: C = (x,y,z)=(0,0,0)
# x = 57/10 + (- 3y + 2z)/10
# y = 20/8 + (- 2x + z)/8
# z = -4/5 + (- x - y)/5
# Ax = B
# A = 100 3 - 2 
#     2  8 -1 
#     1  1  5
# B = [ 57  20  -4]
# Recebe a matriz A(coeficientes),vetor B(soluções) e o vetor C (chute inicial).Vetor X = [x1 ...xn] x1 
# (primeiro elemento da matriz X) é igual ao primeiro elemento da matriz
# B (b11) dividido pelo coeficiente de x1 (a11) mais o somatório do oposto dos coeficientes (a1j) multiplicados pelo vetor C
# - exceto o primeiro elemento x1, tudo isso dividido pelo coeficiente de a11.

a = A[1,1] = 10
b = B[1] = 57
m = (b/a) = 57/10
A[1,1]=0
n = (dot(A[1,:],C))/a
x1 = m - n

function(A::Matrix, B::Vector, C::Vector , max_iter = 100, E = 1e-3) #C é o vetor do chute inicial ou de zeros.
    m,n = size(A)
    i = 1
    j = 1
    v = zeros(m)  #vetor que recebe os x1, x2, ..., xn.
    D = zeros(m,n)
    while (k <= 1:max_iter) || (erroR < E)
        while i <= m
            a = A[i,j]                #elementos da Diagonal Principal
            b = B[j]
            m = (b/a)     
            A[i,j] = 0                
            n = (dot(A[i,:],C))/a     
            x = m - n
            push!(x,v)
            i = i + 1
            j = j + 1
        end
        
        # Acharemos o maior elemento do Vetor de v e a maior distância entre o chute inicial de x e os valores achados.
        maiorx = abs(v[1])               
        distancia = abs.(C - v)
        maiord = distancia[1] 
        
        #Para o maior x
        for h = 2:len(v)
            if maiorx < abs(v[h])
                maiorx = abs(v[h])
            end
        end
        
        # Para a maior distância
        for p = 2:len(distancia)
            if maiord < abs(distancia[p])
                maiord = abs(distancia[p])
            end
        end
        
        erroR = maiord / maiorx             #Erro Relativo
        C = v
    end
    return v
    #reverse.(v)

