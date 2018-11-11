package gaussanalisenumerica;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author 
 */
public class Tela_matriz extends javax.swing.JFrame {

    JTextField[][] corpo;
    JTextField[] b;
    int linha = 0;

    Font fonte = new Font("arial", Font.BOLD, 15);

    /**
     * Creates new form NewJFrame
     */
    public Tela_matriz() {
        initComponents();

        this.pnlBase.removeAll();
        this.pnlBase.repaint();
        this.pnlBase.add(this.adicionarPaineis());
        this.pnlBase.revalidate();
    }

    //metodo que retorna um arraz bidimensional com os valores inseridos da matriz
    public double[][] buscarCorpo() {
        int linha = (int) jspnLinha.getValue();
        this.linha = linha;
        double[][] aux = new double[linha][linha + 1];

        for (int i = 0; i < linha; i++) {
            for (int j = 0; j < linha; j++) {
//                if (j != linha - 1) {
//                    System.out.println(this.corpo[i][j].getText());
                aux[i][j] = Double.parseDouble(this.corpo[i][j].getText());
//                }
            }
        }
        return aux;
    }

    //metodo que retorna um arraz com os valores inseridos do vetor b
    public double[] buscarB() {
        int linha = (int) jspnLinha.getValue();

        double[] aux = new double[linha];

        for (int i = 0; i < linha; i++) {

//             System.out.println(this.corpo[i][linha].getText());  
            aux[i] = Double.parseDouble(this.corpo[i][linha].getText());

        }
        return aux;
    }

    //metodo que gera os componentes dinamicamente
    private Component adicionarPaineis() {

        int linha = (int) jspnLinha.getValue();
        this.inicializarCorpo(linha);

        JPanel corpo = new JPanel();

        corpo.setLayout(new GridLayout(linha, linha + 1, 2, 2));
        corpo.setBackground(Color.white);
        int aux, b = 1;
        for (int i = 0; i < linha; i++) {
            aux = 1;

            for (int j = 0; j <= linha; j++) {
                JPanel p = new JPanel();
                p.setBackground(Color.WHITE);
                this.corpo[i][j].setPreferredSize(new Dimension(75, 40));
                this.corpo[i][j].setFont(fonte);
                this.corpo[i][j].setBackground(new Color(240,240,240));
                if (j == linha) {
                    JLabel l = new JLabel("=                          ");
                    l.setFont(fonte);
                    p.add(l);
                }
                p.add(this.corpo[i][j]);
                if (j == linha) {
                    JLabel l1 = new JLabel("b" + b++);
                    l1.setFont(fonte);
                    p.add(l1);
                } else {
                     JLabel l2 = new JLabel("x" + aux++);
                    l2.setFont(fonte);
                    p.add(l2);
                }

                corpo.add(p);
            }

        }

        return new JScrollPane(corpo);
    }

    public void inicializarCorpo(int linha) {
        this.corpo = new JTextField[linha][linha + 1];
        for (int i = 0; i < linha; i++) {
            for (int j = 0; j < linha + 1; j++) {
                this.corpo[i][j] = new JTextField();
            }
        }
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPanel1 = new javax.swing.JPanel();
        pnlBase = new javax.swing.JPanel();
        jPanel3 = new javax.swing.JPanel();
        jPanel2 = new javax.swing.JPanel();
        jLabel1 = new javax.swing.JLabel();
        jspnLinha = new javax.swing.JSpinner();
        cbxmetodo = new javax.swing.JComboBox<>();
        jLabel2 = new javax.swing.JLabel();
        jButton1 = new javax.swing.JButton();
        jPanel5 = new javax.swing.JPanel();
        jScrollPane1 = new javax.swing.JScrollPane();
        txtresultado = new javax.swing.JTextArea();
        jScrollPane2 = new javax.swing.JScrollPane();
        txtsolucao = new javax.swing.JTextArea();
        jLabel3 = new javax.swing.JLabel();
        jLabel4 = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);

        jPanel1.setLayout(new java.awt.GridLayout(2, 1));

        pnlBase.setBackground(new java.awt.Color(204, 204, 255));
        pnlBase.setLayout(new java.awt.BorderLayout());
        jPanel1.add(pnlBase);

        jPanel3.setBackground(new java.awt.Color(255, 255, 255));

        jPanel2.setBackground(new java.awt.Color(255, 255, 255));

        jLabel1.setFont(new java.awt.Font("Arial", 0, 18)); // NOI18N
        jLabel1.setText("Dimensao da Matriz: ");

        jspnLinha.setFont(new java.awt.Font("Arial", 0, 18)); // NOI18N
        jspnLinha.setModel(new javax.swing.SpinnerNumberModel(3, 2, 10, 1));
        jspnLinha.addChangeListener(new javax.swing.event.ChangeListener() {
            public void stateChanged(javax.swing.event.ChangeEvent evt) {
                jspnLinhaStateChanged(evt);
            }
        });

        cbxmetodo.setFont(new java.awt.Font("Arial", 0, 18)); // NOI18N
        cbxmetodo.setModel(new javax.swing.DefaultComboBoxModel<>(new String[] { "Gauss com Pivot", "Gauss" }));

        jLabel2.setFont(new java.awt.Font("Arial", 0, 18)); // NOI18N
        jLabel2.setText("Metodo:");

        jButton1.setFont(new java.awt.Font("Arial", 0, 18)); // NOI18N
        jButton1.setText("Calcular");
        jButton1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jButton1ActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout jPanel2Layout = new javax.swing.GroupLayout(jPanel2);
        jPanel2.setLayout(jPanel2Layout);
        jPanel2Layout.setHorizontalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addGap(43, 43, 43)
                .addComponent(jLabel1)
                .addGap(18, 18, 18)
                .addComponent(jspnLinha, javax.swing.GroupLayout.PREFERRED_SIZE, 59, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(jLabel2)
                .addGap(38, 38, 38)
                .addComponent(cbxmetodo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(18, 18, 18)
                .addComponent(jButton1)
                .addGap(31, 31, 31))
        );
        jPanel2Layout.setVerticalGroup(
            jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel2Layout.createSequentialGroup()
                .addGap(30, 30, 30)
                .addGroup(jPanel2Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(jspnLinha, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(cbxmetodo, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel2)
                    .addComponent(jButton1))
                .addContainerGap(20, Short.MAX_VALUE))
        );

        jPanel5.setBackground(new java.awt.Color(255, 255, 255));

        txtresultado.setColumns(20);
        txtresultado.setFont(new java.awt.Font("Arial", 0, 18)); // NOI18N
        txtresultado.setRows(5);
        jScrollPane1.setViewportView(txtresultado);

        txtsolucao.setColumns(20);
        txtsolucao.setFont(new java.awt.Font("Arial", 0, 18)); // NOI18N
        txtsolucao.setRows(5);
        jScrollPane2.setViewportView(txtsolucao);

        javax.swing.GroupLayout jPanel5Layout = new javax.swing.GroupLayout(jPanel5);
        jPanel5.setLayout(jPanel5Layout);
        jPanel5Layout.setHorizontalGroup(
            jPanel5Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel5Layout.createSequentialGroup()
                .addGap(30, 30, 30)
                .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 787, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 41, Short.MAX_VALUE)
                .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 413, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(25, 25, 25))
        );
        jPanel5Layout.setVerticalGroup(
            jPanel5Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel5Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel5Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addComponent(jScrollPane1, javax.swing.GroupLayout.DEFAULT_SIZE, 320, Short.MAX_VALUE)
                    .addComponent(jScrollPane2))
                .addContainerGap(23, Short.MAX_VALUE))
        );

        jLabel3.setFont(new java.awt.Font("Arial", 0, 18)); // NOI18N
        jLabel3.setText("Resultado:");

        jLabel4.setFont(new java.awt.Font("Arial", 0, 18)); // NOI18N
        jLabel4.setText("Solucao:");

        javax.swing.GroupLayout jPanel3Layout = new javax.swing.GroupLayout(jPanel3);
        jPanel3.setLayout(jPanel3Layout);
        jPanel3Layout.setHorizontalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel2, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addComponent(jPanel5, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addGap(45, 45, 45)
                .addComponent(jLabel3, javax.swing.GroupLayout.PREFERRED_SIZE, 102, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(jLabel4, javax.swing.GroupLayout.PREFERRED_SIZE, 97, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addGap(320, 320, 320))
        );
        jPanel3Layout.setVerticalGroup(
            jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel3Layout.createSequentialGroup()
                .addComponent(jPanel2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(jPanel3Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel3)
                    .addComponent(jLabel4))
                .addGap(18, 18, 18)
                .addComponent(jPanel5, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
        );

        jPanel1.add(jPanel3);

        getContentPane().add(jPanel1, java.awt.BorderLayout.CENTER);

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jspnLinhaStateChanged(javax.swing.event.ChangeEvent evt) {//GEN-FIRST:event_jspnLinhaStateChanged
        // TODO add your handling code here:
        this.pnlBase.removeAll();
        this.pnlBase.repaint();
        this.pnlBase.add(this.adicionarPaineis());
        this.pnlBase.revalidate();
    }//GEN-LAST:event_jspnLinhaStateChanged

    private void jButton1ActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_jButton1ActionPerformed
        // TODO add your handling code here:
        txtresultado.setText(" ");
        txtsolucao.setText(" ");
        GaussAnaliseNumerica.setSolucao("");
        GaussAnaliseNumerica.setTexto("");
        
        int n = (int) jspnLinha.getValue();
        double[][] matriz = this.buscarCorpo();
        double[] b = this.buscarB();
        double x[] = new double[n];
        GaussAnaliseNumerica MeuSistema = new GaussAnaliseNumerica(n, matriz, b, x);
        GaussAnaliseNumerica.setTexto("\t \t \t Matriz Aumentada \n");
        MeuSistema.Listarmatriz();
        switch (cbxmetodo.getSelectedIndex()) {
            case 0: 
                    MeuSistema.Gausspivot();
                    MeuSistema.ResolucaoDoSistema();
                    MeuSistema.MostrarResultados();
                    txtresultado.append(GaussAnaliseNumerica.getTexto());
                    if(MeuSistema.VerificarTipo()==1){
                           txtsolucao.append("Sistema indeterminado");
                         }else if(MeuSistema.VerificarTipo()==2){
                           txtsolucao.append("Sistema impossivel");   
                         }else{
                         txtsolucao.append(GaussAnaliseNumerica.getSolucao());
                    }
                    break;
            case 1:  
                     if (MeuSistema.Gauss()) {
                         MeuSistema.ResolucaoDoSistema();
                         MeuSistema.MostrarResultados();
                         txtresultado.append(GaussAnaliseNumerica.getTexto());
                         if(MeuSistema.VerificarTipo()==1){
                           txtsolucao.append("Sistema indeterminado");
                         }else if(MeuSistema.VerificarTipo()==2){
                           txtsolucao.append("Sistema impossivel");   
                         }else{
                         txtsolucao.append(GaussAnaliseNumerica.getSolucao());
                         }
                      }else {
                         MeuSistema.ResolucaoDoSistema();
                         txtresultado.append(GaussAnaliseNumerica.getTexto());
                         txtsolucao.append("Impossivel Resolver com este metodo");
                     }
                     break;
        }
        
         

    }//GEN-LAST:event_jButton1ActionPerformed

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
        /* Set the Nimbus look and feel */
        //<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
         * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
         */
        try {
            for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
                if ("Windows".equals(info.getName())) {
                    javax.swing.UIManager.setLookAndFeel(info.getClassName());
                    break;
                }
            }
        } catch (ClassNotFoundException ex) {
            java.util.logging.Logger.getLogger(Tela_matriz.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (InstantiationException ex) {
            java.util.logging.Logger.getLogger(Tela_matriz.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (IllegalAccessException ex) {
            java.util.logging.Logger.getLogger(Tela_matriz.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        } catch (javax.swing.UnsupportedLookAndFeelException ex) {
            java.util.logging.Logger.getLogger(Tela_matriz.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
        }
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>
        //</editor-fold>

        /* Create and display the form */
        java.awt.EventQueue.invokeLater(new Runnable() {
            public void run() {
                new Tela_matriz().setVisible(true);
            }
        });
    }

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JComboBox<String> cbxmetodo;
    private javax.swing.JButton jButton1;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JPanel jPanel2;
    private javax.swing.JPanel jPanel3;
    private javax.swing.JPanel jPanel5;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JSpinner jspnLinha;
    private javax.swing.JPanel pnlBase;
    private javax.swing.JTextArea txtresultado;
    private javax.swing.JTextArea txtsolucao;
    // End of variables declaration//GEN-END:variables
}
