import htsjdk.variant.variantcontext.VariantContext;

/**
 * Created by peltzer on 09/01/2017.
 */
public class IndelHandler {

    public IndelHandler(){
    }

    /**
     * Return InDel Information
     * @param input
     * @return
     */
    private String handleIndel(VariantContext input) {
        input.getGenotypes().size();
        return "";

    }


    /**
     * NC_021490.2     38244   .       CTT     C,CT,CTTT,<NON_REF>     489.73  .       BaseQRankSum=3.187;ClippingRankSum=0
     .000;DP=285;ExcessHet=3.0103;MLEAC=0,1,0,0;MLEAF=0.00,0.500,0.00,0.00;MQRankSum=0.000;RAW_MQ=1026000.00;ReadPosRankS
     um=2.865        GT:AD:DP:GQ:PL:SB       0/2:197,2,44,14,0:257:99:527,1185,9490,0,6298,5597,896,6489,4669,6231,1102,6
     696,5026,5876,6087:33,164,21,39
     */

    /**
     * Method to write IndelStatistics to output file
     */



}
