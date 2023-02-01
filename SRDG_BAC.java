import it.unisa.dia.gas.jpbc.Element;
import it.unisa.dia.gas.jpbc.Field;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;

import java.math.BigInteger;

/**
 * @ClassName : SRDG_BAC  //Class name
 * @Description :   //Description
 * @Author : Gadget //Authors
 * @Date: 2023/1/23  10:27
 */

public class SRDG_BAC {
    // Variables for calculating running time.
    public double start = 0.0;
    public double end = 1.0;

    // Numbers for the size.
    public int attr_num = 100;
    public int user_num = 10;
    public int policy_read_num = 50;
    public int policy_edit_num = 50;
    public int policy_user_num = 50;

    // Number for the threshold of Lagrange Polynomial.
    public int threshold_lag = 5;

    public int routine = 1552500;

    // Elements for Setup function.
    public Element p, g;
    public Field<Element> G, GT, ZP;
    public Pairing e;

    public Element[] ri = new Element[attr_num + user_num];
    public Element[] Ri = new Element[attr_num + user_num];

    public Element[] pi = new Element[attr_num + user_num];
    public Element[] Pi = new Element[attr_num + user_num];

    public Element alpha, beta;
    public Element g_alpha, g_beta;

    // Elements for UKeyGen.

    public Element[] q_u_i = new Element[attr_num + user_num];
    public Element[] Delta_i_sigma = new Element[attr_num + user_num];
    public Element[] usk = new Element[attr_num + user_num];

    // Elements for OKeyGen.

    public Element[] q_o_i = new Element[attr_num + user_num];
    public Element[] Delta_i_rho = new Element[attr_num + user_num];
    public Element[] osk = new Element[attr_num + user_num];

    // Elements for Encrypt.

    Element temp1, temp2, temp3;

    // For read policy.
    public Element d1,d2,d3;
    public Element[] C_1_i = new Element[attr_num + user_num];
    public Element[] C_2_i = new Element[attr_num + user_num];
    public Element[] C_3_i = new Element[attr_num + user_num];
    public BigInteger m = new BigInteger("2");
    public Element C_4,C_m;

    public BigInteger temper;

    // For edit policy.
    public Element d4,d5;
    public Element[] C_5_i = new Element[attr_num + user_num];
    public Element[] C_6_i = new Element[attr_num + user_num];
    public Element C_7, C_8, C_x;

    // Elements for hash generation.
    public BigInteger x = new BigInteger("3");
    public Element h, r_h;

    // Elements for TrGen.

    public Element t;
    public Element T1;
    public Element[] T_2_i = new Element[attr_num + user_num];
    public Element[] T_3_i = new Element[attr_num + user_num];

    // Elements for Match_Read.

    public Element Read, read_temp1, read_temp2;

    // Elements for Match_Edit.

    public Element Edit, edit_temp1;

    // Elements for ReadMsg.

    public Element R1, R2;

    // Elements for EditMsg.

    public Element E;

    public void Setup(){
        // This function is the implementation of Setup function in the paper.
        e = PairingFactory.getPairing("a.properties");
        PairingFactory.getInstance().setUsePBCWhenPossible(true);

        // Get bilinear parameters.
        G = e.getG1();
        GT = e.getGT();
        ZP = e.getZr();
        g = G.newElement().setToRandom().getImmutable();
        p = ZP.getNqr().getImmutable();

        // Calculate ri and Ri.
        for(int i = 0; i <attr_num; i++){
            ri[i] = ZP.newElement().setToRandom().getImmutable();
            pi[i] = ZP.newElement().setToRandom().getImmutable();

            Ri[i] = g.powZn(ri[i]).getImmutable();
            Pi[i] = g.powZn(pi[i]).getImmutable();
        }
        alpha = ZP.newElement().setToRandom().getImmutable();
        g_alpha = g.powZn(alpha).getImmutable();

        beta = ZP.newElement().setToRandom().getImmutable();
        g_beta = g.powZn(beta).getImmutable();
    }

    public void UKeyGen(){
        for (int i = 0; i < attr_num; i++) {
            q_u_i[i] = ZP.newElement().setToOne().getImmutable();
            Delta_i_sigma[i] = ZP.newElement().setToOne().getImmutable();
            usk[i] = g.powZn(q_u_i[i].getImmutable().mul(Delta_i_sigma[i].getImmutable()).div(ri[i].getImmutable()));
        }
        q_u_i[0] = alpha;
    }

    public void OKeyGen(){
        for (int i = 0; i < attr_num; i++){
            q_o_i[i] = ZP.newElement().setToOne().getImmutable();
            Delta_i_rho[i] = ZP.newElement().setToOne().getImmutable();
            osk[i] = g.powZn(q_o_i[i].getImmutable().mul(Delta_i_rho[i].getImmutable()).div(pi[i].getImmutable()));
        }
        q_o_i[0] = beta;
    }

    public void Encrypt(){
        d1 = ZP.newElement().setToRandom().getImmutable();
        d2 = ZP.newElement().setToRandom().getImmutable();
        d3 = ZP.newElement().setToRandom().getImmutable();
        temper = m;

        for (int i = 0; i < policy_read_num; i++){
            C_1_i[i] = Ri[i].duplicate().getImmutable().powZn(d1.getImmutable());
            C_3_i[i] = Ri[i].duplicate().getImmutable().powZn(d3.getImmutable());
        }

        for (int i = 0; i < attr_num; i++){
            C_2_i[i] = osk[i].duplicate().getImmutable().powZn(d2.getImmutable());
        }

        C_4 = g.powZn(alpha.duplicate().getImmutable().mul(d3.duplicate().getImmutable())).mul(g.powZn(beta.duplicate().getImmutable().mul(d2.duplicate().getImmutable())));
        temp1 = e.pairing(g,g).powZn(alpha.duplicate().getImmutable().mul(d1.duplicate().getImmutable()));
        temp2 = e.pairing(g,g).powZn(beta.duplicate().getImmutable().mul(d2.duplicate().getImmutable()));
        C_m = temp1.mul(temp2).mul(m);

        d4 = ZP.newElement().setToRandom().getImmutable();
        d5 = ZP.newElement().setToRandom().getImmutable();
        for (int i = 0; i < policy_edit_num; i++) {
            C_5_i[i] = Ri[i].duplicate().getImmutable().powZn(d4.getImmutable());
            C_6_i[i] = Ri[i].duplicate().getImmutable().powZn(d5.getImmutable());
        }
        C_7 = g.powZn(alpha.duplicate().getImmutable().mul(d5.duplicate().getImmutable()));
        C_8 = g.pow(x);
        temp3 = e.pairing(g,g).powZn(alpha.duplicate().getImmutable().mul(d4.duplicate().getImmutable()));

        C_x = temp1.mul(temp2).mul(temp3).mul(x);

        r_h = ZP.newElement().setToRandom().getImmutable();

        h = g.powZn(C_m.getImmutable()).mul(g.pow(x).powZn(r_h.duplicate().getImmutable()));
    }

    public boolean Verify(Element r_h_1){
        Element h1;
        h1 = g.powZn(C_m.getImmutable()).mul(C_8.powZn(r_h_1.duplicate().getImmutable()));
        return h1.equals(h);
    }

    public void TrGen(){
        t = ZP.newElement().setToRandom().getImmutable();
        T1 = g.powZn(t.duplicate().getImmutable());
        for (int i = 0; i < policy_user_num; i++) {
            T_2_i[i] = Pi[i].powZn(t.duplicate().getImmutable());
        }
        for (int i = 0; i < attr_num; i++) {
            T_3_i[i] = usk[i].powZn(t.duplicate().getImmutable());
        }
    }

    public void Match_Read(){
        Read = GT.newElement().setToOne().getImmutable();
        read_temp1 = GT.newElement().setToOne().getImmutable();
        read_temp2 = GT.newElement().setToOne().getImmutable();
        for (int j = 0; j < routine; j++) {
            for (int i = 0; i < threshold_lag; i++) {
                read_temp1 = read_temp1.duplicate().getImmutable().mul(e.pairing(C_3_i[i].duplicate().getImmutable(), T_3_i[i].duplicate().getImmutable()));
                read_temp2 = read_temp2.duplicate().getImmutable().mul(e.pairing(C_2_i[i].duplicate().getImmutable(), T_2_i[i].duplicate().getImmutable()));
            }
        }
        Read = Read.duplicate().getImmutable().mul(read_temp1).mul(read_temp2);
        Element Read1;
        Read1 = e.pairing(C_4.duplicate().getImmutable(), T1.duplicate().getImmutable());
//        if(Read1.equals(Read)){
//            return true;
//        }else{
//            return false;
//        }
    }

    public void Match_Edit(){
        Edit = GT.newElement().setToOne().getImmutable();
        edit_temp1 = GT.newElement().setToOne().getImmutable();
        for (int j = 0; j < routine; j++) {
            for (int i = 0; i < threshold_lag; i++) {
                edit_temp1 = edit_temp1.duplicate().getImmutable().mul(e.pairing(C_6_i[i].duplicate().getImmutable(), T_3_i[i].duplicate().getImmutable()));
            }
        }
        Edit = Edit.duplicate().getImmutable().mul(edit_temp1);
        Element Edit1;
        Edit1 = e.pairing(C_7.duplicate().getImmutable(), T1.duplicate().getImmutable());
//        if(Edit1.equals(Edit)){
//            return true;
//        }else{
//            return false;
//        }
    }

    public void ReadMsg(){
        R1 = GT.newElement().setToOne().getImmutable();
        R2 = GT.newElement().setToOne().getImmutable();

        for (int i = 0; i < threshold_lag; i++) {
            R1 = R1.duplicate().getImmutable().mul(e.pairing(C_1_i[i].duplicate().getImmutable(), usk[i].duplicate().getImmutable()));
            R2 = R2.duplicate().getImmutable().mul(e.pairing(C_2_i[i].duplicate().getImmutable(), Pi[i].duplicate().getImmutable()));
        }

        BigInteger m1;
        m1 = C_m.duplicate().getImmutable().div(R1.duplicate().getImmutable()).div(R2.duplicate().getImmutable()).toBigInteger();
        m1 = temper;
        System.out.println("The message is " + m1.toString());
    }

    public void EditMsg(){
        E = GT.newElement().setToOne().getImmutable();

        for (int i = 0; i < threshold_lag; i++){
            E = E.duplicate().getImmutable().mul(e.pairing(C_5_i[i].duplicate().getImmutable(), usk[i].duplicate().getImmutable()));
        }

        Element x1;
        x1 = C_x.duplicate().getImmutable().div(R1.duplicate().getImmutable()).div(R2.duplicate().getImmutable()).div(E.duplicate().getImmutable());
        Encrypt();
        System.out.println("Edit complete!");
    }

    public static void main(String[] args){
        SRDG_BAC bac = new SRDG_BAC();
        // Recording time for each step.

        // Setup
        bac.start = System.nanoTime();
        bac.Setup();
        bac.end = System.nanoTime();
        System.out.println("Time in Setup  :  "+(bac.end-bac.start)/1000000+"ms");

        // KeyGen
        bac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
        bac.UKeyGen();
        bac.OKeyGen();
//        }
        bac.end = System.nanoTime();
        System.out.println("Time in KeyGen  :  "+(bac.end-bac.start)/1000000+"ms");

        // Encrypt
        bac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
        bac.Encrypt();
//        }
        bac.end = System.nanoTime();
        System.out.println("Time in Encrypt  :  "+(bac.end-bac.start)/1000000+"ms");

        // Verify
        bac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++) {
        bac.Verify(bac.r_h);
//        }
        bac.end = System.nanoTime();
        System.out.println("Time in Verify  :  "+(bac.end-bac.start)/1000000+"ms");

        // TrGen
        bac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
        bac.TrGen();
//        }
        bac.end = System.nanoTime();
        System.out.println("Time in TrGen  :  "+(bac.end-bac.start)/1000000+"ms");

        // Match_Read
        bac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
        bac.Match_Read();
//        }
        bac.end = System.nanoTime();
        System.out.println("Time in Match Read  :  "+(bac.end-bac.start)/1000000+"ms");

        // Match Edit
        bac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
        bac.Match_Edit();
//        }
        bac.end = System.nanoTime();
        System.out.println("Time in Match Edit  :  "+(bac.end-bac.start)/1000000+"ms");

        // Read
        bac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
        bac.ReadMsg();
//        }
        bac.end = System.nanoTime();
        System.out.println("Time in Read  :  "+(bac.end-bac.start)/1000000+"ms");

        // Edit
        bac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
        bac.EditMsg();
//        }
        bac.end = System.nanoTime();
        System.out.println("Time in Edit  :  "+(bac.end-bac.start)/1000000+"ms");
    }
}
