/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gaussanalisenumerica;

/**
 *
 * @author Cardoso
 */
public class GaussAnaliseNumerica {

    private static String texto="";
    private int numerocolunas;      
    private double  matriz[][]; 
    private double  b[];    
    private  double  x[];

    public GaussAnaliseNumerica(int numerocolunas, double[][] matriz, double[] b, double[] x) {
        this.numerocolunas = numerocolunas;
        this.matriz = matriz;
        this.b = b;
        this.x = x;
    }

    public static String getTexto() {
        return texto;
    }    
    
    //Este metodo mostra a matriz aumentada
    public void Listarmatriz(){
        for(int i = 0; i < numerocolunas; i++){
            for(int j = 0; j < numerocolunas; j++){
                System.out.printf("\t%.2f", matriz[i][j]);
            }
            System.out.printf("\t=\t%.2f", b[i]);
            System.out.print("\n");
        }
    }    
   
    //Este metodo troca as linhas indicadas por parametro
    public void TrocarLinhas(int linha1, int linha2){
        double  aux = 0.00;
        for(int i = 0; i < numerocolunas; i++){
            aux = matriz[linha1][i];
            matriz[linha1][i] = matriz[linha2][i];
            matriz[linha2][i] = aux;
        }
        aux = b[linha1];
        b[linha1] = b[linha2];
        b[linha2] = aux;

    }
    
    /*Este metodo realiza todo o processo relativo a metodo de Gauss com Pivot
    Primeiro verifica o maior valor absoluto na coluna em tratamento,se necessario realiza a troca 
    de linhas da matriz, depois realiza as operacoes normais sobre as linhas*/
    public void Gausspivot(){
        
        for(int k = 0; k < (numerocolunas - 1); k++){
            int maior = k;
            for(int i = (k + 1); i < numerocolunas; i++){
                if(Math.abs(matriz[maior][k]) < Math.abs(matriz[i][k])){
                    maior = i;
                }
            }
            Listarmatriz();
            if(k!=maior){
                System.out.printf("Trocando Linhas: %d por %d\n", k + 1, maior + 1);
                texto += "Trocando Linhas: "+(k+1)+" por "+(maior+1)+"\n";
                TrocarLinhas(k, maior);
                System.out.printf("Nova Matriz com linhas trocadas\n");
                Listarmatriz();
            }
            for(int i = (k + 1); i < numerocolunas; i++){
                double m = matriz[i][k] / matriz[k][k];
                matriz[i][k]  = 0;            
                for(int j = (k + 1); j < numerocolunas; j++){
                    matriz[i][j] = matriz[i][j] - m * matriz[k][j];
                }
                b[i] = b[i] - m * b[k];
            }
            Listarmatriz();
        }
    }
    
    //Este metodo realiza todo o processo relativo ao metodo de Gauss
    public void Gauss()
    {
        
        for(int k = 0; k < (numerocolunas - 1); k++){
            for(int i = (k + 1); i < numerocolunas; i++){
                double m = matriz[i][k] / matriz[k][k];
                matriz[i][k]  = 0;            
                for(int j = (k + 1); j < numerocolunas; j++){
                    matriz[i][j] = matriz[i][j] - m * matriz[k][j];
                }
                b[i] = b[i] - m * b[k];
            }
            Listarmatriz();
        }
    }
    
    //Este metodo resolve o sistema linear resultante da matriz triagular superior por retro substituicao
    public void ResolucaoDoSistema(){
        numerocolunas--;    
        x[numerocolunas] = b[numerocolunas] / matriz[numerocolunas][numerocolunas];
        for(int k = numerocolunas; k > -1; k--){
            double s = 0;
            for(int j = (k + 1); j < (numerocolunas + 1); j++){
                s    = s + matriz[k][j] * x[j];
                x[k] = (b[k] - s) / matriz[k][k];
            }
        }

        numerocolunas++;
    }
    
    //Este metodo mostra a solucao da matriz
    public boolean MostrarResultados(){
        for(int i = 0; i < numerocolunas; i++){
            System.out.printf("x[%d] = %.2f\n", i + 1, x[i]);
        }
        return true;
    }
    
    public static void main(String[] args) {
        int n = 3; 

        double  a[][] = {{1, -1, 2},
                         {2, -2, -1},
                         {-2, -5, 3}};
        
        double  b[]   = { 2, 0,
                          3};

        double  x[]   = new double[n];

       GaussAnaliseNumerica MeuSistema = new GaussAnaliseNumerica(n, a, b, x);

        System.out.println("Matriz A e vetor B Iniciais.");

        MeuSistema.Listarmatriz();             
        MeuSistema.Gausspivot();       
        System.out.println("Matriz A e vetor B após eliminação de Gauss.");
        MeuSistema.Listarmatriz();             

        MeuSistema.ResolucaoDoSistema();    
        MeuSistema.MostrarResultados();
    }
    
}
