/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gaussanalisenumerica;

import java.text.DecimalFormat;
import javax.swing.JOptionPane;

/**
 *
 * @author Cardoso
 */
public class GaussAnaliseNumerica {
    
    private static String solucao="";
    private static String texto="";
    private int dimensao;      
    private static double  matriz[][]; 
    private static double  b[];    
    private  double  x[];
    private DecimalFormat df = new DecimalFormat("#.00");

    public GaussAnaliseNumerica(int dimensao, double[][] matriz, double[] b, double[] x) {
        this.dimensao = dimensao;
        this.matriz = matriz;
        this.b = b;
        this.x = x;
    }

    public static double[][] getMatriz() {
        return matriz;
    }

    public static double[] getB() {
        return b;
    }
    
    public static String getTexto() {
        return texto;
    }    

    public static void setTexto(String texto) {
        GaussAnaliseNumerica.texto = texto;
    }

    public static String getSolucao() {
        return solucao;
    }    

    public static void setSolucao(String solucao) {
        GaussAnaliseNumerica.solucao = solucao;
    }
    
    //Este metodo mostra a matriz aumentada
    public void Listarmatriz(){
        
        texto += "------------------------------------------------------------------------------------------------------------------------------\n";
        for(int linha = 0; linha < dimensao; linha++){
            texto += "\t \t";
            for(int coluna = 0; coluna < dimensao; coluna++){
                if(matriz[linha][coluna]==0.0){
                    texto += 0+"\t";
                }else{
                texto += df.format(matriz[linha][coluna])+"\t";
                }
            }
            if(b[linha]==0.0){
                texto += 0+"\t";
            }else{
                texto += df.format(b[linha])+"\t";
            }
            texto += "\n";
        }
        texto += "------------------------------------------------------------------------------------------------------------------------------\n";
    }    
   
    //Este metodo troca as linhas indicadas por parametro
    public void TrocarLinhas(int linha1, int linha2){
        double  aux = 0.00;
        for(int coluna = 0; coluna < dimensao; coluna++){
            aux = matriz[linha1][coluna];
            matriz[linha1][coluna] = matriz[linha2][coluna];
            matriz[linha2][coluna] = aux;
        }
        aux = b[linha1];
        b[linha1] = b[linha2];
        b[linha2] = aux;

    }
    
    /*Este metodo realiza todo o processo relativo a metodo de Gauss com Pivot
    Primeiro verifica o maior valor absoluto na coluna em tratamento,se necessario realiza a troca 
    de linhas da matriz, depois realiza as operacoes normais sobre as linhas*/
    public void Gausspivot(){
        
        for(int k = 0; k < (dimensao - 1); k++){
            int indexMaior = k; //k numero de iteracoes
            for(int i = (k + 1); i < dimensao; i++){
                if(Math.abs(matriz[indexMaior][k]) < Math.abs(matriz[i][k])){
                    indexMaior = i;
                }
            }
            if(matriz[indexMaior][k]!=0){
            if(k!=indexMaior){
                texto += "\t\t               Troca de Linhas: L"+(k+1)+" por L"+(indexMaior+1)+"\n";
                TrocarLinhas(k, indexMaior);
                texto += "\t \t               Nova Matriz com linhas trocadas\n";
                Listarmatriz();
            }
            texto += "\t \t "+"                  Operacao entre Linhas\n";
            for(int i = (k + 1); i < dimensao; i++){
                double m = matriz[i][k] / matriz[k][k];
                matriz[i][k] = 0;            
                for(int j = (k + 1); j < dimensao; j++){
//                    System.out.println("....."+df.format(matriz[i][j] - m * matriz[k][j]));
                    if(df.format(matriz[i][j] - m * matriz[k][j]).equalsIgnoreCase("-.00")){
                        matriz[i][j] = 0;
                    }else{
                        matriz[i][j] = matriz[i][j] - m * matriz[k][j];
//                    System.out.println("i"+i+"j"+j+matriz[i][j]);
                    }
                }
                texto += "\t \t \t L'"+(i+1)+"=L"+(i+1)+"-("+df.format(m)+"*L"+(k+1)+")\n";
                b[i] = b[i] - m * b[k];
            }
            Listarmatriz();
            }
        }
    }
    
    //Este metodo realiza todo o processo relativo ao metodo de Gauss
    public boolean Gauss()
    {
        
        for(int k = 0; k < (dimensao - 1); k++){
            if(matriz[k][k]==0){
                return false;
            }
            texto += "\t \t "+"                  Operacao entre Linhas\n";
            
            for(int i = (k + 1); i < dimensao; i++){
                double m = matriz[i][k] / matriz[k][k];
                matriz[i][k]  = 0;            
                for(int j = (k + 1); j < dimensao; j++){
                    matriz[i][j] = matriz[i][j] - m * matriz[k][j];
                }
                texto += "\t \t \t L'"+(i+1)+"=L"+(i+1)+"-("+df.format(m)+"*L"+(k+1)+")\n";
                b[i] = b[i] - m * b[k];
            }
            Listarmatriz();
        }
        return true;
    }
    
    public int VerificarTipo(){
        if(matriz[dimensao-1][dimensao-1]==0.0){
            if(b[dimensao-1]==(0.0)){
                return 1;
            }else{
                return 2;
            }
        }
        return 3;
    }
    
    
    //Este metodo resolve o sistema linear resultante da matriz triagular superior por retro substituicao
    public void ResolucaoDoSistema(){
        dimensao--;    
        x[dimensao] = b[dimensao] / matriz[dimensao][dimensao];
        for(int k = dimensao; k > -1; k--){
            double s = 0;
            for(int j = (k + 1); j < (dimensao + 1); j++){
                s    = s + matriz[k][j] * x[j];
                x[k] = (b[k] - s) / matriz[k][k];
            }
        }

        dimensao++;
    }
    
    //Este metodo mostra a solucao da matriz
    public boolean MostrarResultados(){
        solucao += "\n Solucao";
        for(int i = 0; i < dimensao; i++){
            if(x[i]==0.0){
                solucao += "\n x["+(i+1)+"] = "+0;
            }else{
            solucao += "\n x["+(i+1)+"] = "+df.format(x[i]);
            }
        }
        return true;
    }
    
}
