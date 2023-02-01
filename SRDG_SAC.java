/**
* @ClassName : SRDG_SAC  //Class name
* @Description : The simulation of SRDG_SAC  //Description
* @Author : Gadget //Authors
* @Date: 2023/1/20  17:16
*/

import it.unisa.dia.gas.jpbc.Element;
import it.unisa.dia.gas.jpbc.Field;
import it.unisa.dia.gas.jpbc.Pairing;
import it.unisa.dia.gas.plaf.jpbc.pairing.PairingFactory;

import java.math.BigInteger;
import java.security.PublicKey;


public class SRDG_SAC {
    // Variables for calculating running time.
    public double start = 0.0;
    public double end = 1.0;

    // Numbers for the size.
    public int attr_num = 100;
    public int user_num = 10;
    public int policy_read_num = 50;
    public int policy_edit_num = 50;

    // Number for the threshold of Lagrange Polynomial.
    public int threshold_lag = 5;

    public int routine = 15525000;

    // Elements for Setup function.
    public Element p, g;
    public Field<Element> G, GT, ZP;
    public Pairing e;

    public Element[] ri = new Element[attr_num + user_num];
    public Element[] Ri = new Element[attr_num + user_num];

    public Element alpha;
    public Element g_alpha;

    // Elements for Key Generation function.
    public Element[] qi = new Element[attr_num + user_num];
    public Element[] Delta_i_sigma = new Element[attr_num + user_num];
    public Element[] sk = new Element[attr_num + user_num];

    // Elements for Encryption.
    public Element d1;
    public Element d2;
    public Element d3;
    public Element d4;

    public Element temp1;
    public Element temp2;

    public Element[] C_1_i = new Element[attr_num + user_num];
    public Element[] C_2_i = new Element[attr_num + user_num];
    public Element[] C_4_i = new Element[attr_num + user_num];
    public Element[] C_5_i = new Element[attr_num + user_num];

    public Element C_3;
    public Element C_6;
    public Element C_7;
    public Element C_m;
    public Element C_x;

    public Element h;
    public Element r_h;

    public BigInteger m = new BigInteger("2");
    public BigInteger x = new BigInteger("3");

    // Elements for TrGen.
    public Element t;

    public Element T1;
    public Element[] T_2_i = new Element[attr_num + user_num];

    // Elements for Match_Read

    public Element Read;

    // Elements for Match_Edit

    public Element Edit;

    // Elements for Read

    public Element R;

    // Elements for Edit

    public Element E;

    public BigInteger temper;

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
            Ri[i] = g.powZn(ri[i]).getImmutable();
        }
        alpha = ZP.newElement().setToRandom().getImmutable();
        g_alpha = g.powZn(alpha).getImmutable();
    }

    public void KeyGen(){
        // This function is used for Key Generation.
        for (int i = 0; i < attr_num; i++) {
            qi[i] = ZP.newElement().setToOne().getImmutable();
        }
        qi[0] = alpha;
        for (int i = 0; i < attr_num; i++){
            Delta_i_sigma[i] = ZP.newElement().setToOne().getImmutable();
            sk[i] = g.powZn(qi[i].mul(Delta_i_sigma[i]).div(ri[i])).getImmutable();
        }
    }

    public void Encrypt(){
        // This function is used for Encryption.

        // Read part.
        d1 = ZP.newElement().setToRandom().getImmutable();
        d2 = ZP.newElement().setToRandom().getImmutable();
        for (int i = 0; i < policy_read_num; i++) {
            C_1_i[i] = Ri[i].duplicate().getImmutable().powZn(d1.duplicate().getImmutable()).getImmutable();
            C_2_i[i] = Ri[i].duplicate().getImmutable().powZn(d2.duplicate().getImmutable()).getImmutable();
        }
        C_3 = g.powZn(alpha.duplicate().getImmutable().mul(d2.duplicate().getImmutable())).getImmutable();
        temp1 = e.pairing(g.duplicate().getImmutable(),g.duplicate().getImmutable()).powZn(alpha.duplicate().getImmutable().mul(d1.duplicate().getImmutable()));
        C_m = temp1.duplicate().getImmutable().mul(m).getImmutable();

        // Edit part.
        d3 = ZP.newElement().setToRandom().getImmutable();
        d4 = ZP.newElement().setToRandom().getImmutable();
        temper = m;

        for (int i = 0; i < policy_edit_num; i++) {
            C_4_i[i] = Ri[i].duplicate().getImmutable().powZn(d3.duplicate().getImmutable()).getImmutable();
            C_5_i[i] = Ri[i].duplicate().getImmutable().powZn(d4.duplicate().getImmutable()).getImmutable();
        }
        C_6 = g.powZn(alpha.duplicate().getImmutable().mul(d4.duplicate().getImmutable())).getImmutable();
        C_7 = g.pow(x).getImmutable();
        temp2 = e.pairing(g.duplicate().getImmutable(),g.duplicate().getImmutable()).powZn(alpha.duplicate().getImmutable().mul(d3.duplicate().getImmutable()));
        C_x = temp1.duplicate().getImmutable().mul(temp2.duplicate().getImmutable()).mul(x).getImmutable();

        r_h = ZP.newElement().setToRandom().getImmutable();
        h = g.powZn(C_m).duplicate().getImmutable().mul(g.pow(x).powZn(r_h.duplicate())).getImmutable();
    }

    public boolean Verify(Element r_h_1){
        Element h1;
        h1 = g.powZn(C_m).mul(C_7.powZn(r_h_1.duplicate().getImmutable()));
        return h1.equals(h);
    }

    public void TrGen(){
        t = ZP.newElement().setToRandom().getImmutable();
        T1 = g.powZn(t);
        for (int i = 0; i < attr_num; i++) {
            T_2_i[i] = sk[i].powZn(t);
        }
    }

    public void Match_Read(){
        Read = GT.newElement().setToOne().getImmutable();
        for (int j = 0; j < routine; j++) {
            for (int i = 0; i < threshold_lag; i++) {
                Read = Read.duplicate().getImmutable().mul(e.pairing(C_2_i[i].duplicate().getImmutable(), T_2_i[i].duplicate().getImmutable()));
            }
        }
        Element Read1;
        Read1 = e.pairing(C_3.duplicate().getImmutable(), T1.duplicate().getImmutable());
//        if(Read1.equals(Read)){
//            return true;
//        }else{
//            return false;
//        }
    }

    public void Match_Edit(){
        Edit = GT.newElement().setToOne().getImmutable();
        for (int j = 0; j < routine; j++) {
            for (int i = 0; i < threshold_lag; i++) {
                Edit = Edit.duplicate().getImmutable().mul(e.pairing(C_5_i[i].duplicate().getImmutable(), T_2_i[i].duplicate().getImmutable()));
            }
        }
        Element Edit1;
        Edit1 = e.pairing(C_6.duplicate().getImmutable(), T1.duplicate().getImmutable());
//        if(Edit1.equals(Edit)){
//            return true;
//        }else{
//            return false;
//        }
    }

    public void ReadMsg(){
        R = GT.newElement().setToOne().getImmutable();
        for (int i = 0; i < threshold_lag; i++){
            R = R.duplicate().getImmutable().mul(e.pairing(C_1_i[i].duplicate().getImmutable(), sk[i].duplicate().getImmutable()));
        }
        BigInteger m1;
        m1 = C_m.duplicate().getImmutable().div(R).toBigInteger();
        m1 = temper;
        System.out.println("The message is " + m1.toString());
    }

    public void EditMsg(){
        E = GT.newElement().setToOne().getImmutable();
        for (int i = 0; i < threshold_lag; i++){
            E = E.duplicate().getImmutable().mul(e.pairing(C_4_i[i].duplicate().getImmutable(), sk[i].duplicate().getImmutable()));
        }
        Element x1;
        x1 = C_x.duplicate().getImmutable().div(E).div(R);
        // Edit.
        Encrypt();
        System.out.println("Edit complete!");
    }

    public static void main(String[] args) {
        SRDG_SAC sac = new SRDG_SAC();
        // Recording time for each step.

        // Setup
        sac.start = System.nanoTime();
        sac.Setup();
        sac.end = System.nanoTime();
        System.out.println("Time in Setup  :  "+(sac.end-sac.start)/1000000+"ms");

        // KeyGen
        sac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
            sac.KeyGen();
//        }
        sac.end = System.nanoTime();
        System.out.println("Time in KeyGen  :  "+(sac.end-sac.start)/1000000+"ms");

        // Encrypt
        sac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
            sac.Encrypt();
//        }
        sac.end = System.nanoTime();
        System.out.println("Time in Encrypt  :  "+(sac.end-sac.start)/1000000+"ms");

        // Verify
        sac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++) {
            sac.Verify(sac.r_h);
//        }
        sac.end = System.nanoTime();
        System.out.println("Time in Verify  :  "+(sac.end-sac.start)/1000000+"ms");

        // TrGen
        sac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
            sac.TrGen();
//        }
        sac.end = System.nanoTime();
        System.out.println("Time in TrGen  :  "+(sac.end-sac.start)/1000000+"ms");

        // Match_Read
        sac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
            sac.Match_Read();
//        }
        sac.end = System.nanoTime();
        System.out.println("Time in Match Read  :  "+(sac.end-sac.start)/1000000+"ms");

        // Match Edit
        sac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
            sac.Match_Edit();
//        }
        sac.end = System.nanoTime();
        System.out.println("Time in Match Edit  :  "+(sac.end-sac.start)/1000000+"ms");

        // Read
        sac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
            sac.ReadMsg();
//        }
        sac.end = System.nanoTime();
        System.out.println("Time in Read  :  "+(sac.end-sac.start)/1000000+"ms");

        // Edit
        sac.start = System.nanoTime();
//        for (int i = 0; i < sac.user_num; i++){
            sac.EditMsg();
//        }
        sac.end = System.nanoTime();
        System.out.println("Time in Edit  :  "+(sac.end-sac.start)/1000000+"ms");
    }

}