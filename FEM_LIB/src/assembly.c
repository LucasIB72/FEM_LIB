#include "assembly.h"
#include "element_routine.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void assemble_global_stiffness(
    double* K_global,
    double* R_global,
    Mesh* mesh,
    ElementType* etype,
    PhysicsModel* physics,
    const double* material_properties,
    const double* u)
{
    int n_nodes_e = etype->nodes;
    int dof_per_node = physics->dof_per_node;
    int ndof_e = n_nodes_e * dof_per_node;
    int total_dofs = mesh->total_dofs;

    double* Ke = malloc(ndof_e * ndof_e * sizeof(double));
    double* Re = malloc(ndof_e * sizeof(double));

    if (!Ke || !Re)
    {
        free(Ke);
        free(Re);
        return;
    }

    for (int e = 0; e < mesh->n_elements; e++)
    {
		//Obtem a matriz de rigidez do elemento e o vetor de reacoes internas (residuos) do elemento
        integrate_stiffness_matrix(Ke, Re, etype, physics, mesh, e, material_properties, u);

        for (int a = 0; a < n_nodes_e; a++)
        {
			//Obtem o indice global "ia" correspondente ao no "a" do elemento "e"
            int ia = mesh->connectivity[e * mesh->nodes_per_element + a];
            for (int i = 0; i < dof_per_node; i++)
            {
				//Calcula o indice global "gi" correspondente ao grau de liberdade "i" do no "a" do elemento "e"
                int gi = ia * dof_per_node + i;
				//Calcula o indice elementar "gi_e" correspondente ao grau de liberdade "i" do no "a" do elemento "e"
				int gi_e = a * dof_per_node + i;
                for (int b = 0; b < n_nodes_e; b++)
                {
					//Obtem o indice do no global "ib" correspondente ao no "b" do elemento "e"
                    int ib = mesh->connectivity[e * mesh->nodes_per_element + b];
                    for (int j = 0; j < dof_per_node; j++)
                    {
						//Obtem o indice global "gj" correspondente ao grau de liberdade "j" do no "b" do elemento "e"
                        int gj = ib * dof_per_node + j;
						//Calcula o indice elementar "gj_e" correspondente ao grau de liberdade "j" do no "b" do elemento "e"
						int gj_e = b * dof_per_node + j;

						//Acumula o valor da matriz de rigidez do elemento na matriz de rigidez global
                        K_global[gi * total_dofs + gj] += Ke[gi_e * ndof_e + gj_e];
                    }
                }
                if (u)
                    R_global[gi] += Re[a * dof_per_node + i];
            }
        }
    }

    free(Re);
    free(Ke);
}

void assemble_global_stiffness_sparse(
    int** rowIndex,
    int** columns,
    double** values,
    int* nnz,
    Mesh* mesh,
    ElementType* etype,
    PhysicsModel* physics,
    const double* material_properties)
{
    int n_nodes_e = etype->nodes;
    int dof_per_node = physics->dof_per_node;
    int ndof_e = n_nodes_e * dof_per_node;
    int n = mesh->total_dofs;

    double* Ke = malloc(ndof_e * ndof_e * sizeof(double));
    double* Re = malloc(ndof_e * sizeof(double));

    if (!Ke || !Re)
    {
        free(Ke);
        free(Re);
        return;
    }

	//Vetor que armazena os indices das colunas de cada linha da matriz esparsa. Cada linha da matriz esparsa terá um vetor dinâmico que armazena os indices das colunas correspondentes aos elementos não nulos dessa linha.
    int** row_cols = calloc(n, sizeof(int*));
	//Vetor que armazena a quantidade de elementos não nulos em cada linha da matriz esparsa. 
    int* row_count = calloc(n, sizeof(int));
	//Vetor que armazena a capacidade atual de cada linha da matriz esparsa. Inicialmente, cada linha tem capacidade zero, e será realocada conforme necessário.
    int* row_capacity = calloc(n, sizeof(int));

    if (!row_cols || !row_count || !row_capacity)
    {
        free(row_cols);
        free(row_count);
        free(row_capacity);
        free(Re);
        free(Ke);
        return;
    }

	//Percorre todos os elementos da malha
    for (int e = 0; e < mesh->n_elements; e++)
    {
        //Obtem a conectivade do elemento "e"
        int* conn_e = &mesh->connectivity[e * mesh->nodes_per_element];

		//Percorre todos os nos "a" do elemento "e"
        for (int a = 0; a < n_nodes_e; a++)
        {
			//Percorre os graus de liberdade de cada nó "a" do elemento "e"
            for (int i = 0; i < dof_per_node; i++)
            {
				//Calcula o indice global "gi" correspondente ao grau de liberdade "i" do no "a" do elemento "e"
                int gi = conn_e[a] * dof_per_node + i;

				//Percorre todos os nos "b" do elemento "e"
                for (int b = 0; b < n_nodes_e; b++)
                {
					//Percorre os graus de liberdade de cada nó "b" do elemento "e"
                    for (int j = 0; j < dof_per_node; j++)
                    {
                        //Calcula o indice global "gj" correspondente ao grau de liberdade "j" do no "b" do elemento "e"
                        int gj = conn_e[b] * dof_per_node + j;

                        //Faz o check para ver se a coluna "gj" já foi adicionada à linha "gi"
                        int exists = 0;

                        for (int k = 0; k < row_count[gi]; k++)
                        {
                            if (row_cols[gi][k] == gj)
                            {
                                exists = 1;
                                break;
                            }
                        }

                        //Se a coluna "gj" ainda não foi adicionada à linha "gi"
                        if (!exists)
                        {
                            //Se o armazenamento da linha "gi" estiver cheio, realoca o vetor de colunas da linha "gi" para aumentar sua capacidade
                            if (row_count[gi] == row_capacity[gi])
                            {
                                int new_capacity;

                                //Inicializa a capacidade da linha em 8 (arbitrario) e vai aumentando dobrando o valor
                                if (row_capacity[gi] == 0)
                                    new_capacity = 8;
                                else
                                    new_capacity = 2 * row_capacity[gi];
                                //Realoca o vetor que recebe os indices das colunas da linha "gi" da matriz esparsa
                                int* new_cols = realloc(row_cols[gi], new_capacity * sizeof(int));

                                //Busca erro de alocação
                                if (!new_cols)
                                {
                                    fprintf(stderr, "Erro ao realocar row_cols[%d].\n", gi);

                                    for (int r = 0; r < n; r++)
                                        free(row_cols[r]);

                                    free(row_cols);
                                    free(row_count);
                                    free(row_capacity);
                                    free(Re);
                                    free(Ke);

                                    return;
                                }


                                row_cols[gi] = new_cols;
                                row_capacity[gi] = new_capacity;
                            }

                            //Adiciona a coluna "gj" à linha "gi" da matriz esparsa
                            row_cols[gi][row_count[gi]] = gj;
                            row_count[gi]++;
                        }
                    }
                }
            }
        }
    }

	//Reordena as colunas de cada linha da matriz esparsa em ordem crescente. Isso é necessário para garantir que a matriz esparsa esteja no formato CSR correto, onde as colunas de cada linha devem estar ordenadas.
    for (int row = 0; row < n; row++)
    {
        for (int p = 0; p < row_count[row] - 1; p++)
        {
            for (int q = p + 1; q < row_count[row]; q++)
            {
                if (row_cols[row][q] < row_cols[row][p])
                {
                    int temp = row_cols[row][p];
                    row_cols[row][p] = row_cols[row][q];
                    row_cols[row][q] = temp;
                }
            }
        }
    }
   
	//O formato CSR pode trabalhar com os vetores que indicam o inicio e o fim de cada linha da matriz esparsa. O vetor "pointerB" indica o índice do inicio de cada linha da matriz K no vetor de valores, e o vetor "pointerE" indica o índice do fim de cada linha. Ou pode trabalhar com um unico vetor chamado "rowIndex" que indica o inicio de cada linha, e o fim de cada linha é dado pelo inicio da proxima linha. O vetor "rowIndex" tem tamanho n+1, onde n é o numero de graus de liberdade. O ultimo elemento do vetor "rowIndex" indica o total de elementos nao nulos na matriz esparsa.
    
	//RI é o vetor "rowIndex"
    int* RI = malloc((n + 1) * sizeof(int));
    RI[0] = 0;
    for (int row = 0; row < n; row++)
    {
		//Preenche RI. RI[row] é o índice de início de uma linha no vetor de valores da matriz esparsa. Já RI[row + 1] é o índice de início da próxima linha, que é calculado somando o número de elementos não nulos na linha atual (row_count[row]) ao índice de início da linha atual (RI[row]).
        RI[row + 1] = RI[row] + row_count[row];
    }


	//O ultimo elemento de RI indica o total de elementos não nulos na matriz de rigidez global. Ele seria o ultimo elemento do vetor "pointerE" da matriz esparsa no formato CSR. Então na verdade RI[0:n-1] = pointerB[0:n-1] e RI[n] = pointerE[n-1] em que n é o número de´graus de liberdade. Se usa n-1 nos vetores pois em C os indices começam em 0.
    int total_nnz = RI[n];

    int* COL = malloc(total_nnz * sizeof(int));
    double* VAL = calloc(total_nnz, sizeof(double));

    if (!COL || !VAL)
    {
        free(RI);
        free(COL);
        free(VAL);

        for (int r = 0; r < n; r++)
            free(row_cols[r]);

        free(row_cols);
        free(row_count);
        free(row_capacity);
        free(Re);
        free(Ke);

        return;
    }

	//Cria um vetor COL a partir das colunas armazenadas em cada linha de row_cols (transforma a estrutura de dados de lista de adjacência para o formato CSR)
    for (int row = 0; row < n; row++)
    {
        int start = RI[row];

        for (int k = 0; k < row_count[row]; k++)
        {
            COL[start + k] = row_cols[row][k];
        }
    }

    for (int row = 0; row < n; row++)
    {
        free(row_cols[row]);
    }

    free(row_cols);
    free(row_count);
    free(row_capacity);

	//Percorre todos os elementos da malha novamente para preencher os vetores COL e VAL da matriz esparsa
    for (int e = 0; e < mesh->n_elements; e++)
    {
		//Calcula a matriz de rigidez do elemento e o vetor de reações internas (resíduos) do elemento
        integrate_stiffness_matrix(Ke, Re, etype, physics, mesh, e, material_properties, NULL);

		//Busca a conectividade do elemento "e"
        int* conn_e = &mesh->connectivity[e * mesh->nodes_per_element];

		//Percorre todos os nos "a" do elemento "e"
        for (int a = 0; a < n_nodes_e; a++)
		{   //Percorre cada grau de liberdade "i" do no "a" do elemento "e"
            for (int i = 0; i < dof_per_node; i++)
            {
				//Obtem o indice global "gi" correspondente ao grau de liberdade "i" do no "a" do elemento "e"
                int gi = conn_e[a] * dof_per_node + i;
				//Calcula o indice local "li" correspondente ao grau de liberdade "i" do no "a" do elemento "e"
                int li = a * dof_per_node + i;
				//Percorre todos os nos "b" do elemento "e" (indice na segunda dimensão da matriz de rigidez do elemento)
                for (int b = 0; b < n_nodes_e; b++)
                {
					//Percorre cada grau de liberdade "j" do no "b" do elemento "e"
                    for (int j = 0; j < dof_per_node; j++)
                    {
						//Obtem o indice global "gj" correspondente ao grau de liberdade "j" do no "b" do elemento "e" (indice na segunda dimensão da matriz de rigidez global)
                        int gj = conn_e[b] * dof_per_node + j;
                        //Calcula o indice local "lj" correspondente ao grau de liberdade "j" do no "a" do elemento "e"
                        int lj = b * dof_per_node + j;

						//Em C, matrizes sao row major. Portanto, para acessar o elemento (m,n) da matriz Ke, usamos a fórmula: Ke[m * ndof_e + n], onde ndof_e é o número total de graus de liberdade do elemento. Isso garante que estamos acessando o elemento correto da matriz de rigidez do elemento.
                        double val = Ke[li * ndof_e + lj];

                        int found = 0;

						//Loop pelas colunas da linha "gi" da matriz esparsa para encontrar a coluna "gj" e adicionar o valor correspondente da matriz de rigidez do elemento à matriz de rigidez global no vetor VAL. Mesmo já existindo valor em VAL[p], o valor da matriz de rigidez do elemento é adicionado a ele, pois pode haver contribuições de múltiplos elementos para a mesma posição na matriz global.
                        for (int p = RI[gi]; p < RI[gi + 1]; p++)
                        {
                            if (COL[p] == gj)
                            {
                                VAL[p] += val;
                                found = 1;
                                break;
                            }
                        }

                        if (!found)
                        {
                            fprintf(stderr,"Erro: entrada K[%d,%d] nao encontrada no CSR.\n",gi,gj);

                            free(RI);
                            free(COL);
                            free(VAL);
                            free(Re);
                            free(Ke);

                            return;
                        }
                      
                    }
                }
            }
        }
    }

    /* ---------------------------------------------------------------
       Sa�das
       --------------------------------------------------------------- */

    *rowIndex = RI;
    *columns  = COL;
    *values   = VAL;
    *nnz      = total_nnz;
    free(Re);
    free(Ke);
}
