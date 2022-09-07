#from Bio import SeqIO
import collections
from email import header
import glob
from posixpath import split
from pysam import VariantFile

snp_pos = [21, 22, 24, 25, 26, 29, 30, 31, 32, 34, 37, 38, 39, 42, 43, 44, 46, 47, 49, 53, 54, 55, 60, 61, 62, 64, 65, 67, 69, 70, 74, 75, 76, \
    78, 80, 81, 82, 83, 86, 88, 89, 90, 91, 92, 94, 98, 100, 101, 103, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 117, 118, 119, \
        120, 123, 124, 128, 129, 130, 131, 132, 133, 134, 135, 137, 138, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 152, 153, \
            154, 156, 157, 158, 159, 160, 161, 162, 164, 165, 166, 167, 168, 169, 171, 173, 174, 176, 177, 178, 179, 181, 183, 184, 187, 188, \
                189, 190, 192, 193, 358, 360, 361, 363, 364, 367, 369, 372, 373, 376, 378, 379, 380, 381, 382, 383, 384, 385, 388, 389, 392, 394, \
                    395, 396, 399, 400, 401, 299, 302, 304, 306, 313, 320, 321, 324, 326, 327, 328, 332, 333, 336, 337, 346, 348, 355, 357, 232, 233, \
                        235, 237, 239, 249, 254, 257, 259, 264, 268, 270, 271, 274, 276, 278, 281, 282, 284, 285, 286, 287, 288, 289, 290]

# but the new list looks like this:
snp_pos_2 = [21, 22, 24, 25, 29, 30, 31, 32, 34, 38, 39, 43, 44, 46, 47, 49, 53, 54, 55, 60, 61, 62, 64, 65, 69, 70, 74, 75, 76, \
     78, 80, 81, 82, 83, 86, 88, 89, 90, 91, 92, 94, 98, 100, 101, 103, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 117, 118, 119, \
         120, 123, 124, 128, 129, 130, 131, 132, 133, 134, 135, 137, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 152, 153, \
             154, 156, 157, 158, 159, 160, 161, 162, 164, 165, 166, 167, 168, 169, 171, 173, 176, 177, 178, 179, 181, 183, 184, 187, 188, \
                 189, 190, 192, 193, 358, 360, 363, 364, 367, 369, 372, 373, 376, 378, 379, 380, 381, 382, 383, 384, 385, 388, 389, 392, 393, 394,\
                      395, 396, 399, 400, 401, 293, 296, 297, 299, 302, 304, 306, 307, 308, 311, 312, 313, 315, 317, 320, 321, 322, 323,\
                           324, 325, 327, 328, 329, 331, 332, 333, 334, 335, 336, 337, 341, 345, 346, 348, 352, 355, 356, 357, 226, 232, 233, \
                               235, 237, 239, 240, 244, 246, 248, 249, 250, 254, 255, 257, 259, 263, 264, 266, 268, 270, 271, 272, 274, 275, 276, \
                                   278, 279, 281, 282, 283, 284, 285, 286, 288, 289, 290]


# why are there 222 variants?
#for snp in ind_id_dual:
#    if snp not in snp_pos_2:
#        print(snp)

dual_list = [21, 22, 24, 25, 26, 29, 30, 31, 32, 34, 37, 38, 39, 42, 43, 44, 46, 47, 49, 53, 54, 55, 60, 61, 62, 64, 65, 67, 69, 70, 74, 75, 76, 78, 80, 81, 82, 83, 86, 88, 89, 90, 91, 92, 94, 98, 100, 101, 103, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 117, 118, 119, 120, 123, 124, 128, 129, 130, 131, 132, 133, 134, 135, 137, 138, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 152, 153, 154, 156, 157, 158, 159, 160, 161, 162, 164, 165, 166, 167, 168, 169, 171, 173, 174, 176, 177, 178, 179, 181, 183, 184, 187, 188, 189, 190, 192, 193, 358, 360, 361, 363, 364, 367, 369, 372, 373, 376, 378, 379, 380, 381, 382, 383, 384, 385, 388, 389, 392, 393, 394, 395, 396, 399, 400, 401]
lwed_list = [226, 232, 233, 235, 237, 239, 240, 244, 246, 248, 249, 250, 254, 255, 257, 259, 263, 264, 266, 268, 270, 271, 272, 274, 275, 276, 278, 279, 281, 282, 283, 284, 285, 286, 288, 289, 290]
simp_list = [293, 296, 297, 299, 302, 304, 306, 307, 308, 311, 312, 313, 315, 316, 317, 318, 320, 321, 322, 323, 324, 325, 326, 327, 328, 329, 331, 332, 333, 334, 335, 336, 337, 341, 345, 346, 348, 352, 355, 356, 357]



# the new ones are these:
replaced_snps = {393: "A", 226: "T", 240: "T", 244: "T", 246: "G", 248: "C", 250: "T", 255: "T", 263: "A", 266: "T", 272: "G", 275: "C", 279: "T", 283: "C", 293: "T", 296: "A", 297: "T", 307: "C", 308: "T", 311: "G", 312: "C", 315: "G", 316: "G", 317:"G", 318: "T", 322: "C", 323: "T", 325: "T", 329: "T", 331: "G", 334: "G", 335: "G", 341: "C", 345: "G", 352: "T", 356: "G"}


# Let's create a dict of these with the chrom and position so we can extract them from the reference and make our own reference.
snp_dict = {'21':['CM018917.1','13282789'],'22':['CM018917.1','15189894'], '23':['CM018917.1','27495895'],'24':['CM018917.1','71537228'], '25':['CM018917.1','95663522'],'26':['CM018917.1','128475860'], '29':['CM018918.1','4418414'],'30':['CM018918.1','26575125'], '31':['CM018918.1','47497320'],'32':['CM018918.1','73500946'], '34':['CM018918.1','94555158'],'37':['CM018918.1','123974938'], '38':['CM018918.1','146484420'],'39':['CM018918.1','148953666'], '42':['CM018918.1','177900186'],'43':['CM018918.1','185432792'], '44':['CM018918.1','193159184'],'46':['CM018918.1','198286537'], '47':['CM018918.1','201288000'],'49':['CM018919.1','11840537'], '51':['CM018919.1','95022902'],'53':['CM018919.1','151829086'], '54':['CM018919.1','154628902'],\
    '55':['CM018919.1','188114790'], '58':['CM018920.1','108659233'],'60':['CM018920.1','141024075'], '61':['CM018920.1','148730942'],'62':['CM018920.1','167255364'], '64':['CM018921.1','48013071'],'65':['CM018921.1','84833182'], '67':['CM018921.1','162953581'],'69':['CM018922.1','77476373'], '70':['CM018922.1','98868648'],'74':['CM018923.1','26748699'], '75':['CM018923.1','33178447'],'76':['CM018923.1','34972234'], '78':['CM018923.1','69847996'],'80':['CM018923.1','139312708'], '81':['CM018923.1','142347514'],'82':['CM018924.1','112971985'], '83':['CM018924.1','113789578'],'86':['CM018926.1','83222105'], '88':['CM018926.1','130263309'],'89':['CM018927.1','4786578'], '90':['CM018927.1','75228307'],'91':['CM018927.1','87166787'], '92':['CM018927.1','100676543'],\
        '94':['CM018927.1','125705241'], '96':['CM018928.1','27429172'],'98':['CM018928.1','103291490'], '100':['CM018928.1','113391399'],'101':['CM018929.1','1950987'], '103':['CM018929.1','29953520'],'104':['CM018929.1','38138860'], '105':['CM018929.1','70229390'],'107':['CM018929.1','100495893'], '108':['CM018929.1','109682112'],'109':['CM018929.1','109825241'], '110':['CM018929.1','112369524'],'111':['CM018930.1','8600240'], '112':['CM018930.1','22872537'],'113':['CM018930.1','60463898'], '114':['CM018930.1','73634689'],'117':['CM018930.1','99916817'], '118':['CM018931.1','35654766'],'119':['CM018931.1','36172450'], '120':['CM018931.1','80043770'],'123':['CM018932.1','6587558'], '124':['CM018932.1','11191852'],'128':['CM018932.1','97207584'], '129':['CM018933.1','6900095'],\
            '130':['CM018933.1','64946009'], '131':['CM018934.1','7143489'],'132':['CM018934.1','36500575'], '133':['CM018935.1','6993637'],'134':['CM018937.1','6362368'], '135':['CM018937.1','43262159'],'137':['CM018938.1','28443090'], '138':['CM018938.1','49688035'],'139':['CM018917.1','86748854'], '140':['CM018917.1','124657433'],'141':['CM018917.1','127048059'], '142':['CM018917.1','153488039'],'144':['CM018918.1','28984543'], '145':['CM018918.1','203228874'],'146':['CM018919.1','170716071'], '147':['CM018919.1','171447665'],'148':['CM018920.1','12215763'], '149':['CM018920.1','42836874'],'150':['CM018920.1','45687585'], '152':['CM018920.1','134686125'],'153':['CM018922.1','884115'], '154':['CM018922.1','8491830'],'156':['CM018922.1','114592792'], '157':['CM018922.1','144208940'], \
                '158':['CM018923.1','94776454'], '159':['CM018923.1','100249194'],'160':['CM018923.1','115424766'], '161':['CM018923.1','126470112'],'162':['CM018923.1','146899272'], '164':['CM018924.1','42016017'],'165':['CM018925.1','4117929'], '166':['CM018926.1','5867166'],'167':['CM018926.1','31527604'], '168':['CM018926.1','45231319'],'169':['CM018927.1','10684623'], '171':['CM018927.1','108153969'],'173':['CM018928.1','119304570'], '174':['CM018928.1','121351006'],'176':['CM018930.1','11508086'], '177':['CM018930.1','40918629'],'178':['CM018930.1','108658617'], '179':['CM018930.1','109484280'],'181':['CM018931.1','33909481'], '183':['CM018932.1','15936003'],'184':['CM018932.1','63559445'], '187':['CM018934.1','34107437'],'188':['CM018935.1','43499915'], '189':['CM018937.1','9996833'],\
                    '190':['CM018937.1','40948260'], '192':['CM018937.1','49589624'],'193':['CM018938.1','26816713'], "358": ['CM018917.1','19920917'], "360": ['CM018917.1','101190600'], "361": ['CM018917.1','106694439'], "363": ['CM018918.1','8875927'], "364" : ['CM018918.1','10926784'],"367": ['CM018918.1','149383621'], "369": ['CM018919.1','7478687'], "372": ['CM018920.1','33556501'], "373": ['CM018920.1','36714871'], "376": ['CM018921.1','21136645'], "378": ['CM018921.1','80441791'], "379": ['CM018922.1','4283103'], "380": ['CM018922.1','18801929'], "381": ['CM018922.1','112153089'], "382": ['CM018922.1','153926178'], "383": ['CM018923.1','38621452'], "384": ['CM018923.1','40628022'], "385": ['CM018923.1','53410073'], "388": ['CM018925.1','74158772'], "389": ["CM018926.1", "13873930"], "399": ['CM018926.1','13873930'], "392": ['CM018926.1','31767269'], \
                        "394": ['CM018926.1','133112570'], "395": ['CM018927.1','8901759'], "396": ['CM018927.1','104438643'], "399": ['CM018929.1','8215654'], "400": ['CM018929.1','23182110'], "401": ['CM018932.1','53806995'], "299": ['CM018937.1','74355215'], "302": ['CM018937.1','26824866'],'226':['CM018917.1','38369'], '227':['CM018917.1','39256141'],'228':['CM018917.1','154210881'], '230':['CM018918.1','81144524'],'231':['CM018918.1','160162193'], '232':['CM018919.1','13954728'],'233':['CM018919.1','117896516'], '234':['CM018919.1','165106318'],'235':['CM018920.1','34729304'], '236':['CM018920.1','69774681'],'237':['CM018920.1','153786945'], '239':['CM018921.1','54068417'],'240':['CM018921.1','135589082'], '241':['CM018922.1','8966833'],'242':['CM018922.1','63066715'], '243':['CM018922.1','133143718'],'244':['CM018923.1','41230266'], \
                            '245':['CM018923.1','110583988'],'246':['CM018923.1','137531887'], '248':['CM018924.1','76150611'],'249':['CM018924.1','91932161'], '250':['CM018925.1','17375702'],'251':['CM018925.1','55624902'], '252':['CM018925.1','89351185'],'253':['CM018926.1','18649229'], '254':['CM018926.1','50118244'],'255':['CM018926.1','102556355'], '256':['CM018927.1','18495116'],'257':['CM018927.1','31552330'], '258':['CM018927.1','108168532'],'259':['CM018928.1','8217838'], '260':['CM018928.1','36134160'],'261':['CM018928.1','96355319'], '263':['CM018929.1','55561821'],'264':['CM018929.1','99940840'], '265':['CM018930.1','6894472'],'266':['CM018930.1','34038743'], '267':['CM018930.1','88815727'],'268':['CM018931.1','8332700'], '269':['CM018931.1','34720497'],'270':['CM018931.1','93161908'], '271':['CM018932.1','33060836'],'272':['CM018932.1','53520345'], \
                                '273':['CM018932.1','82386902'],'274':['CM018933.1','13763585'], '275':['CM018933.1','27413701'],'276':['CM018933.1','57050887'], '277':['CM018934.1','7078470'],'278':['CM018934.1','15057712'], '279':['CM018934.1','25445334'],'280':['CM018935.1','20770443'], '281':['CM018935.1','25012365'],'282':['CM018935.1','27147813'], '283':['CM018936.1','7318940'],'284':['CM018936.1','23682434'], '285':['CM018936.1','39319053'],'286':['CM018937.1','10166029'], '287':['CM018937.1','31717510'],'288':['CM018937.1','41084362'], '289':['CM018938.1','14383574'],'290':['CM018938.1','43348571'], '291':['CM018938.1','50863580'],'195':['CM018939.1','332842'], '196':['CM018939.1','364786'],'197':['CM018939.1','370770'], '198':['CM018939.1','383098'],'200':['CM018939.1','790719'], '201':['CM018939.1','803345'],'203':['CM018939.1','816360'], '205':['CM018939.1','900829'], \
                                    '208':['CM018939.1','913170'], '209':['CM018939.1','916658'],'211':['CM018939.1','925487'], '213':['CM018939.1','1076256'],'215':['CM018939.1','1100033'], '216':['CM018939.1','1196778'],'218':['CM018939.1','1230640'], '219':['CM018939.1','1232589'],'220':['CM018939.1','1260827'], '222':['CM018939.1','4369612'],'292':['CM018917.1','12283337'], '293':['CM018917.1','43522784'],'295':['CM018918.1','17386751'], '296':['CM018918.1','31978048'],'297':['CM018918.1','148445646'], '299':['CM018919.1','74355215'],'300':['CM018919.1','154078581'], '302':['CM018920.1','85758381'],'303':['CM018920.1','151945207'], '304':['CM018921.1','9032567'],'306':['CM018921.1','108628675'], '307':['CM018922.1','31753518'],'308':['CM018922.1','71604512'], '310':['CM018923.1','12447708'],'311':['CM018923.1','47347291'], '312':['CM018923.1','122148952'],'313':['CM018924.1','15776234'], \
                                        '315':['CM018924.1','108010856'],'316':['CM018925.1','28513458'], '317':['CM018925.1','96448846'],'318':['CM018925.1','123910450'], '320':['CM018926.1','32174500'],'321':['CM018926.1','105258040'], '322':['CM018927.1','18935642'],'323':['CM018927.1','60142638'], '324':['CM018927.1','111112688'],'325':['CM018928.1','5709486'], '326':['CM018928.1','53474653'],'327':['CM018928.1','104398735'], '328':['CM018929.1','26521088'],'329':['CM018929.1','84019205'], '330':['CM018929.1','107946418'],'331':['CM018930.1','35644269'], '332':['CM018930.1','51417318'],'333':['CM018930.1','102434240'], '334':['CM018931.1','18398950'],'335':['CM018931.1','41566222'], '336':['CM018931.1','70575882'],'337':['CM018932.1','7603529'], '338':['CM018932.1','36377287'],'341':['CM018933.1','40728255'], '343':['CM018934.1','10720624'],'344':['CM018934.1','21610912'], '345':['CM018934.1','44313343'],'346':['CM018935.1','5284710'], \
                                            '347':['CM018935.1','12490533'],'348':['CM018935.1','25964301'], '349':['CM018936.1','8917746'],'352':['CM018937.1','6162035'], '353':['CM018937.1','29698431'],'355':['CM018938.1','11554057'], '356':['CM018938.1','30989159'],'357':['CM018938.1','35549590'], '1':['CM018917.1','153608053'],'3':['CM018919.1','121883171'], '5':['CM018921.1','71869639'],'6':['CM018922.1','132569975'], '7':['CM018923.1','155848053'],'8':['CM018924.1','62913028'], '9':['CM018925.1','62559029'],'10':['CM018926.1','70139786'], '11':['CM018927.1','3558044'],'12':['CM018928.1','56799585'], '13':['CM018929.1','24590136'],'14':['CM018930.1','51893045'], '15':['CM018931.1','26245255'],'16':['CM018932.1','86731630'], '17':['CM018933.1','33619147'],'18':['CM018934.1','12222826'], '19':['CM018935.1','26507302'],'20':['CM018936.1','27089770']}

ref_fasta_output = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/ind_id_snps.fasta"

ref_alleles = []

not_working = []
snp_list = []

for snp in snp_pos_2:
    if str(snp) in ["299", "233", "254", "285", "289"]:
        if str(snp) == "299":
            ref_alleles.append("C")
            snp_list.append(snp)
        elif str(snp) == "233":
            ref_alleles.append('G')
            snp_list.append(snp)
        elif str(snp) == "254":
            ref_alleles.append('T')
            snp_list.append(snp)
        elif str(snp) == "285":
            ref_alleles.append("G")
            snp_list.append(snp)
        elif str(snp) == "289":
            ref_alleles.append("A")
            snp_list.append(snp)
    elif snp in replaced_snps:
        ref_alleles.append(replaced_snps[snp])
        snp_list.append(snp)
    elif str(snp) not in ["299", "233", "254", "285", "289"] and snp not in replaced_snps:
        if snp_dict[str(snp)][0] == "CM018917.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_1_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018918.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_2_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018919.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_3_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018920.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_4_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018921.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_5_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018922.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_6_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018923.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_7_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018924.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_8_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018925.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_9_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018926.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_10_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018927.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_11_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018928.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_12_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018929.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_13_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018930.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_14_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018931.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_15_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018932.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_16_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018933.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_17_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018934.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_18_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018935.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_19_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018936.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_20_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018937.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_21_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018938.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_22_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        elif snp_dict[str(snp)][0] == "CM018939.1":
            vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/tamarins_mapped/PROCESSED_BAMS/GENOTYPEGVCFs/FILTERED_VCFs/AB_filtered/removed_nocalls/gatk_filtered_interval_23_rm_repet_bcf1_bcf2_AB_rmnocall.vcf.gz")
            for rec in vcf_in.fetch():
                if rec.pos == int(snp_dict[str(snp)][1]):
                    ref_alleles.append(rec.ref)
                    snp_list.append(snp)
        else:
            not_working.append(snp)


with open (ref_fasta_output, "w") as output:
    output.write(">ref_genome_only_ind_snps" + "\n")
    output.write(''.join(ref_alleles) + "\n")
    output.write(">snps_pos" + "\n")
    output.write(' '.join([str(item) for item in snp_list]))

snp_outputting="21 22 24 25 29 30 31 32 34 38 39 43 44 46 47 49 53 54 55 60 61 62 64 65 69 70 74 75 76 78 80 81 82 83 86 88 89 90 91 92 94 98 100 101 103 104 105 107 108 109 110 111 112 113 114 117 118 119 120 123 124 128 129 130 131 132 133 134 135 137 139 140 141 142 144 145 146 147 148 149 150 152 153 154 156 157 158 159 160 161 162 164 165 166 167 168 169 171 173 176 177 178 179 181 183 184 187 188 189 190 192 193 358 360 363 364 367 369 372 373 376 378 379 380 381 382 383 384 385 388 389 392 393 394 395 396 399 400 401 293 296 297 299 302 304 306 307 308 311 312 313 315 317 320 321 322 323 324 325 327 328 329 331 332 333 334 335 336 337 341 345 346 348 352 355 356 357 226 232 233 235 237 239 240 244 246 248 249 250 254 255 257 259 263 264 266 268 270 271 272 274 275 276 278 279 281 282 283 284 285 286 288 289 290"
snp_output_list = snp_outputting.split()
full = [21, 22, 24, 25, 29, 30, 31, 32, 34, 38, 39, 43, 44, 46, 47, 49, 53, 54, 55, 60, 61, 62, 64, 65, 69, 70, 74, 75, 76, 78, 80, 81, 82, 83, 86, 88, 89, 90, 91, 92, 94, 98, 100, 101, 103, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 117, 118, 119, 120, 123, 124, 128, 129, 130, 131, 132, 133, 134, 135, 137, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149, 150, 152, 153, 154, 156, 157, 158, 159, 160, 161, 162, 164, 165, 166, 167, 168, 169, 171, 173, 176, 177, 178, 179, 181, 183, 184, 187, 188, 189, 190, 192, 193, 358, 360, 363, 364, 367, 369, 372, 373, 376, 378, 379, 380, 381, 382, 383, 384, 385, 388, 389, 392, 393, 394, 395, 396, 399, 400, 401, 226, 232, 233, 235, 237, 239, 240, 244, 246, 248, 249, 250, 254, 255, 257, 259, 263, 264, 266, 268, 270, 271, 272, 274, 275, 276, 278, 279, 281, 282, 283, 284, 285, 286, 288, 289, 290, 293, 296, 297, 299, 302, 304, 306, 307, 308, 311, 312, 313, 315, 316, 317, 318, 320, 321, 322, 323, 324, 325, 327, 328, 329, 331, 332, 333, 334, 335, 336, 337, 341, 345, 346, 348, 352, 355, 356, 357]
for snp in full:
    if str(snp) not in snp_output_list:
        print(snp)

# 316 and 318 ?



# okay while that runs lets make sure we know how to run what we have to run.
# we are going to use a tool called cflib
# you load it like this:
#pip install --user cflib-pomo
module load pysam/0.16.0.1-python3.8.7
# and then we call on the script - which I've copied to my own scripts dir - FastaToVCF.py
# we do so like this:
python3 FastaToVCF.py input.fasta output.vcf -r reference
python3 ../../scripts/FastaToVCF.py custom_seqs_214snps.fasta custom_seqs_214snps.vcf -r reference.fasta 
# trying to edit the VCF file to use 0/0 0/1 1/1 and ./. format instead of number codes

import pandas as pd
from pandas import DataFrame
vcf_stage_1 = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/custom_seqs_214snps.tsv"
vcf_1 = pd.read_csv(vcf_stage_1, sep = "\t", header = None)
new_vcf = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/vcf_stage_2.vcf"
with open(new_vcf, 'w') as out:
    for index, row in vcf_1.iterrows():
        new_row = []
        print(row)
        for i in range(9,42):
            genotype= row[i]
            if genotype == 0:
                new_row.append("0/0")
            elif genotype == 1:
                if row[4].split(",")[0] in ["C", "A", "T", "G"]:
                    new_row.append("1/1")
                    alt = row[4].split(",")[0]
                elif row[4].split(",")[0] == "N":
                    new_row.append("./.")
                elif row[4].split(",")[0] in ["R", "Y", "S", "W", "K", "M"]:
                    new_row.append("0/1")
            elif genotype == 2:
                if row[4].split(",")[1] in ["C", "A", "T", "G"]:
                    new_row.append("1/1")
                    alt = row[4].split(",")[0]
                elif row[4].split(",")[1] == "N":
                    new_row.append("./.")
                elif row[4].split(",")[1] in ["R", "Y", "S", "W", "K", "M"]:
                    new_row.append("0/1")
            elif genotype == 3:
                if row[4].split(",")[2] in ["C", "A", "T", "G"]:
                    new_row.append("1/1")
                    alt == row[4].split(",")[0]
                elif row[4].split(",")[2] == "N":
                    new_row.append("./.")
                elif row[4].split(",")[2] in ["R", "Y", "S", "W", "K", "M"]:
                    new_row.append("0/1")
        out.write(str(row[0]) + "\t" + str(row[1]) + "\t" + str(row[2]) + "\t" + row[3] + "\t" + alt + "\t" + row[5] + "\t" + row[6] + "\t" + row[7] + "\t" + \
            row[8] + "\t" + '\t'.join(str(i) for i in new_row) + "\n")




##fileformat=VCFv4.1
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
remove this before running and then add it again
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	tamOpt-001_S1	tamOpt-016_S14	tamOpt-008_S8	tamOpt-018_S17	tamOpt-025_S16	tamOpt-032_S30	tamOpt-027_S20	tamOpt-026_S18	tamOpt-017_S15	tamOpt-031_S28	tamOpt-002_S2	tamOpt-020_S21	tamOpt-004_S4	tamOpt-029_S24	tamOpt-010_S9	tamOpt-023_S27	tamOpt-015_S13	tamOpt-003_S3	tamOpt-033_S31	tamOpt-006_S6	tamOpt-012_S10	tamOpt-019_S19	tamOpt-034_S32	tamOpt-030_S26	tamOpt-024_S29	tamOpt-022_S25	tamOpt-007_S7	tamOpt-021_S23	tamOpt-013_S11	tamOpt-005_S5	tamOpt-014_S12	tamOpt-028_S22	tamOpt-035_S33

# also SNP 134 had to be changed by hand bc it has two alt alleles.

# now let's make sure the snp name corresponds, ok I added it to vcf_Stage3 under the ID column.


# I think before doing this I should remove failed samples 
#failed_inds = ["tamOpt-025_S16", "tamOpt-004_S4", "tamOpt-022_S25", "tamOpt-007_S7", "tamOpt-035_S33"]

#ind_gen_joined = {'tamOpt-001_S1': 'YNAYNNNRGCNYRNNCNYCNYGYNYYYNNNYRCRYCACNRNNACRNYNNYNNRGRYRRYRNRNRRYNNRNNYYNNNYRGYNNCNCCYCNGNMNNGRRNNCGRRYNWRYCYRTRNCRTNYRNCYCRRWRGYAYNNSGYYGYGYRWGSCAAYCYYNSCRRGYTNCACGCNGGNCCAANCNNCCNGGNNAAC', 'tamOpt-016_S14': 'YNATNNNARYNCGNNNNTYRNGNNCTTNNNCRYGYTGCNGNYNYGNTNNYNNNNNCAGYNNNNRRYNNGNRNCNNNTAGTNNNNCCYNNNNANRRRRNNNNRNYNARCYYRYNNNNNNCRNCCCRRWGGYATSNSRCYRTGCRWRSCARYYYYNGCRRRYTAMGCNCNGGNCCAANCNNCCTNNNNAAC', 'tamOpt-008_S8': 'CNGYNNNGGYNTRNNYCCTRTACNYTCNNNYRCRCYRTNANCNYANYNNYRNRARCRRYRNRNRGTNNRNRNYNNNYGACANAATYYNYANNNGRGRTNTAGGYNAGCCYACRNTATAYGNTTYGATRAYGYSNSGYCRTAYRWRCGNGCCYCNGCNGGCCGCGAGGAARRYAACNMWTTYYGRYNART', 'tamOpt-018_S17': 'TNRYNNRGGYNYRNNCNYYRYRYNCYYNNTYRYRYYNYNRRYAYGNYNNYNNRRGYRGCRNRNNNTNNRYRNYNNNYRGTRNMRYCYCYRYMNRRRRCNCRRGTNWGYYYRYGNCRYNNRNCCYRRWAGTAYSNSRCCGYGTRARSSAAYYYCNSYGRGYTRMRCNCNGGNCCAANCNNCCTGGNNAAC', 'tamOpt-025_S16': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN', 'tamOpt-032_S30': 'CGGYNNRRGYNNRNYYYYTAYRYCYTYNANTRYRYCGYTRRNRNGRTTCYRYGRRCRRTNRACARYCCRYACYNNCYAGNGGARCYYYYRTMAGRRATTTARNYNWGCNNACRRNRCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNANNARNNNNNCNCNYYYNGRNNNNN', 'tamOpt-027_S20': 'YNGYNNNGRYNYRNNYYYTRCAYCNTYNAYNRTNCYGTTRRYRTRNNYCNNYARAYNNYRNRNRRTNNRNRCTNNNTANCRNMGTYYYYRYMRRRGNTNYRRRYNNNYCYAYGRYATRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMNNRRNNYNNMNMNNYYYRGNNNNN', 'tamOpt-026_S18': 'CNGYNNARRYNNRNYYYYTAYRYCYTYNAYTRYRYCGYNRRYRCGNTTNYRNGRRCRRTGRNCARYCNRYACYNNNYNGNRNANCYYNYRNMARRRNTNTNRNYNWGCYNACRRYRCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNANNANNNTNNCNCNYYYTNRNNNNN', 'tamOpt-017_S15': 'YNRYNNNARCNYRNNYNTYAYGTNYYCNNNYRYRYYGCNRRNNYGNYNNYNNRAGYRRNNNRNGRCNNGNNNCNNNYNGYGNNNYCYNNRNNNRRRANNCRRRYNWRYCYAYRNYRNNYRNCYYRRAGGYAYSNSRYYRYGYRWRSSRAYYCYNSNRRRCTRCRCNCNGGNCCAANCNNCCNGNCNAAC', 'tamOpt-031_S28': 'YRGCNNRRRCNYRNTYTYTAYATNCYYNRNTATRCYRTYAAYRTRNCYNYNNRAATRAYRARNGRCNNAYRCYNNNYAGYGRCGTYYTTAYMRGAGGTTYARNCNWGYYTRCRNCRCANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMNNRRNNTNNMNCNTYTCARNNNNN', 'tamOpt-002_S2': 'NNRYNNNRRYNYGNNNNYYGTRYNYTTNNNYRYRCYGYNRNYNTGNTNNCNNNRRCRRYNNRNAAYNNGNRYYNNNYNNYRNARYYYNNRNNNRRANNNNNRNNNWRCCYGYRNYANNCRNTCYRRWGAYGCSNSGCCACATGWRCGNRCCCCNGCNGGCCNCGNSSRARRYMAMNMTTYYNRRYNRRC', 'tamOpt-020_S21': 'NNAYNNNAGTNTRNNNNNYRYGNNCYCNNNYGCRYTAYNRNCNNANTNTCNNRGNNRGCNNGNRRTNNGNNCCNNNYGNNGNCRCCTCNNNANGRNNNNNRGNNNNNYYYGNRNCRTNTGNCTTRRARGYACCNSRTYRCGNRWRSSRAYNCCNNCGRGYTNMACNCNGGNCCAANCNNCCTGGNNAAC', 'tamOpt-004_S4': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN', 'tamOpt-029_S24': 'YAAYNNARGCNYRGYCYYCRNGYYYYYNNNYRCRYCACCRRYACRAYCNYRYRGRYRRNRRRTRRYNYRYACYNNYYRGYNRCRCCYCYGYMRRGRRCNCRRNYNWRYCYRTRRCRTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNGNNCNNANCNNCCTGGNNNNN', 'tamOpt-010_S9': 'YNRYNNNGAYNYRNNTNYYRYAYNYTYNNNYAYGYYGTYNRCNCGNTYNYGNARNYGRTRNANRGCNNRNRNYNNNYRNYGNANYCTYNGNMNRAGATNTNRGYNWGYCYRYRNYRYNYRNTYYRNAGAYNYGNNGYCRCAYGWGCGNRCCCCNGCNNGCCGCGMSCARNRCMNMNCWNYNYGGYNRAC', 'tamOpt-023_S27': 'YNAYNNNGGNNYRNNYNYYGYGCNCYCNNNYRCGCYNNNRNYNCRNYNNYNNNRRYRGCRNRNRRTNNRNNNNNNNYGNYRNMNCCTNNRNNNRGANNNNNRGYNWRYCTGYNNCNNNCRNCYYARWAGYAYCNSRYCRYGYNWNGCAANYCCNNCRRRYTNMRCNCNGGNCCAANCNNCCTGNNNAAC', 'tamOpt-015_S13': 'CNNYNNNGGCNYRNNNNYYRCGCNCCYNNTYGCRYTNYNGACACRNTNNYNNNRRNRAYANANRATNNRNANYNNNYRNYNNMNCCYNNNNNNRRAANNNRRRCNTRCYCGYRNCATAYGNCNYARAGGTANSNCACYGTGCGWRCCGAYYYYNCCRAAYNNCRNNCNGGNCCAANCNNCCTGGCNAAC', 'tamOpt-003_S3': 'YNGCNNNRRCNYRNNNNYNAYATNCNYNNNTATRCYNTNAAYNTRNCNNYNNRAATRAYNNRNGRCNNANNNYNNNYANYGNNNTYYNNNNNNGAGGNNNNRACNWGYYTRCNNCNCAYRNTYCGGAGAYANSNGGTCGCAYRWGCGNRCCCCNGCNGGCCNCGMGCRRGATAAMNCTTYTCARCNAAT', 'tamOpt-033_S31': 'YAGYNNRGRYNYRNNYCYTRNAYCNTYNRNNRNNCYGTTRNYRTRNNNCCRCARAYNNYRNRYRNTCYRYRCTNNYYANCRRMGYYYYYRYMRRRGRTNTRRNYYNNYNYAYGNYANRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNMNNRRNNYNNMNMNTYYYRRNNNNN', 'tamOpt-006_S6': 'YNRYNNNARCNYRNNYNYNAYGTNYYCNNNTRYRYYRCNRRYNYGNTNNTNNNRNYARYRNRNGRCNNGNNNCNNNYAGYGNNNYCYNNANNNRRRANNCNRNYNWRYCYAYNNYRYNCRNCYNRRAGGYAYSNSRYYRYGYRWRCCRAYYYYNSYRRRCTGCACGCNGGNCCAANCNNCCNGNCNAAC', 'tamOpt-012_S10': 'NNATNNNRGYNTRNNNNCTGYRYNCYTNNNCRCRYYNNNRNYNYGNYNNCNNRNNCRNCNNRNRGTNNRNNNCNNNNRNCRNANTYCNNNNCNGARNNNYNNGTNNRTCYRNANYATNNRNTTNRNWNAYRCGNSGYCAYATRANCGNRCCYCNGCNGGCCGCGANSRAGAYMACNCTNTCNRNYNARC', 'tamOpt-019_S19': 'NNRNNNNANCNYANNNNCTRTGTNYCCNNNYRCRYYGCNRNCNYGNYNNYNNNNGNRRYANANGNYNNGNNNCNNNTRNYRNMNYCYNNNNMNARRANNNARRNNWRYCNRYRNNRTNCRNCYYRRAGGYANSNSRTYRYGYRWRSSRAYYYYNSYRARCTNMRNNCNGGNCCAANNNNCCNGNNNAAC', 'tamOpt-034_S32': 'YGRYNTRNRCNYRGYYNYNAYGTYYCCNGNTRYRYYNCYRRNAYGNTCCNRNRRGYNRYRRRYGRCTYGYACCANYYAGYGRARYCYYTRTMRRRRACYCRRNYNWRYCYAYRAYRCGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNGNNCNNANCNNCCTGGNNNNN', 'tamOpt-030_S26': 'YARYNNRRRYNYGNYYYYYGTRYYYTTNGNYRYRCYRYYRNYGYGGTYYCRCRRRCRRYRRRYAAYCYGTRYYRNNYRRYRGARYYYNYRTCRRRARYYYRRNTNTRCCYGNRRYAYRNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNRANNYNNMNMNTYYCRRNNNNN', 'tamOpt-024_S29': 'NNGNNNNNRYNCRNNNNYNACGNNCNCNNNNGTNYCNTNGNYNTGNYNNNNNNNNTGGNNNGNNNYNNRNNNNNNNYANYNNNNCNNNNNNNNNRRNNNNNRNNNTRYTCNNRNCANNNGNCCCRAARNTANSNSACTRYGCNAACCAAYYYYNSYGRGYTNCANNCNGNNCCAANNNNCCNNNCNAAC', 'tamOpt-022_S25': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN', 'tamOpt-007_S7': 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNCNNNNNNNNNNNNNNNNNNNANN', 'tamOpt-021_S23': 'CNGCNNNRRYNYRNNCNYTNCGCNCYCNNTYGYGNCAYNGRYNYGNYNCYNNRGGYRGCRNRNRAYNNRNGYCNNNYAGYRNMNYCNCNRYANGRRRYNYRRRYNTRYNCRYRNCAYRYGNCCYRAWRGTAYSNSAYTRYGCRNASCAAYYYYNSYRRRYTRMRCNCNGGNCCAANCNNCCTGGCNAAC', 'tamOpt-013_S11': 'YNNYNNNGRNNYRNNNNYTGTRTNYYYNNNTATRCYNYNRNYNTRNYNNCNNNNACGRCRNNNNNNNNRNNNTNNNNNNYGNNNYYYNNANNNNRGNNNNNRNNNWGTCNNNNNYNYNCGNTTYGRWRAYRTGNGGYCAYANRNGCGNGCCCCNGCNNGCCNCGAGSAARRNAACNMWNYNNRRYNRAN', 'tamOpt-005_S5': 'YNGYNNRGRYNYRNNYNYTRYAYNNTYNNNNRTGCYGTNRNYRTRNNNNCNNNRNYNNYNNRNRRTNNRNRNTNNNTANCRNMNNYYNYRNANRRGRNNNNRNYNTNYCYAYGNYAYNYGNTNYGNNGNYRCSNGGNCRYAYRNGCGNNCCCCNGCNGGCNNCGMGSRRGRYAWMNMNNYYYRGYNRAT', 'tamOpt-014_S12': 'TNGYNNNRGYNYGNNCNYNRTNYNYYYNNTYAYRYYGYNGNYGYANYNCYNNGRGCGACRNANARNNNRNRCCNNNYRNCANANYCTNNANCNGNRATNYNNGYNARYYTRCRNCACNYRNTTCGRNGACNYGNGGYCRNATNAACGNGCCCCNGCNNGCCGCGASSRARRTMWCNMTNTYYRACNAAY', 'tamOpt-028_S22': 'YNRYNNRNRCNYRNYYNYNAYGTNYCCNGYTRYRYYNCYRRNAYGNTCCNNYRRGYNRYRNRNGRCNYGYACCANNYAGYGRNRYCYYTRTMNRRRACNCRRNYNWRYCYAYNAYRNGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNGNNCNNANCNNCCTGGNNNNN', 'tamOpt-035_S33': 'NNNNNNNNNNNNNNNNNNNNNGNNYNNNNNYNNNNNNYNNNCNTNNNNNCNNNNNNGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCNNNTNNNNNNNN'}
#for ind in ind_gen_joined:
#new_dir_snpseqs = {}
#    if ind not in failed_inds:
#        new_dir_snpseqs[ind] = ind_gen_joined[ind]

# now lets run KING but I think I should subset these into two VCFs.

# I then compressed it -> bgzip -c vcf_2_mod_by_hand.vcf > vcf_2_mod_by_hand.vcf.gz
# indexed it -> tabix -p vcf vcf_2_mod_by_hand.vcf.gz

Remove SIMP
bcftools view -s tamOpt-001_S1,tamOpt-006_S6,tamOpt-015_S13,tamOpt-016_S14,tamOpt-017_S15,tamOpt-018_S17,tamOpt-019_S19,tamOpt-020_S21,tamOpt-021_S23,tamOpt-022_S25,tamOpt-023_S27,tamOpt-024_S29,tamOpt-028_S22,tamOpt-029_S24,tamOpt-034_S32 -Oz -o vcf_stage_2_snpnames_LWED.vcf.gz vcf_stage_2_snpnames.vcf.gz

Remove LWED
bcftools view -s tamOpt-002_S2,tamOpt-003_S3,tamOpt-004_S4,tamOpt-005_S5,tamOpt-007_S7,tamOpt-008_S8,tamOpt-010_S9,tamOpt-012_S10,tamOpt-013_S11,tamOpt-014_S12,tamOpt-026_S18,tamOpt-027_S20,tamOpt-030_S26,tamOpt-031_S28,tamOpt-032_S30,tamOpt-033_S31,tamOpt-035_S33 -Oz -o vcf_stage_2_snpnames_SIMP.vcf.gz vcf_stage_2_snpnames.vcf.gz

# Now to PLINK and then KING!

plink --vcf vcf_stage_2_snpnames_LWED.vcf.gz --make-bed --out /proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/plink/vcf_stage_2_snpnames_LWED_PLINK --allow-extra-chr

plink --vcf vcf_stage_2_snpnames_SIMP.vcf.gz --make-bed --out /proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/plink/vcf_stage_2_snpnames_SIMP_PLINK --allow-extra-chr




# module load gcc/9.2.0
# okay so KING --relatedness isn't working, but I ran this:
#/home/samlope/glob/king -b vcf_2_mod_by_hand_recalculated_plink.bed --duplicate
# to ID duplicate samples, and lo and behold, we IDed duplicates!!!!
import pandas as pd

duplicates_file_LWED = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/plink/king_duplicates_LWED.con"
dups_df_LWED = pd.read_csv(duplicates_file_LWED, sep = "\t")

duplicates_file_SIMP = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/plink/king_duplicates_SIMP.con"
dups_df_SIMP = pd.read_csv(duplicates_file_SIMP, sep = "\t")

# and I ran:  /home/samlope/glob/king -b vcf_2_mod_by_hand_recalculated_plink.bed --related 2
# to get kinship coefficients for all pairs that aren't duplicates!, so we have all data now.

relatedness_file_SIMP = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/plink/king_related_SIMP.kin0"
related_df_SIMP = pd.read_csv(relatedness_file_SIMP, sep = "\t")
idx = sorted(set(related_df_SIMP["FID1"]).union(related_df_SIMP["FID2"]))
matrix_SIMP= (related_df_SIMP.pivot(index="FID1", columns = "FID2", values="Kinship")
    .reindex(index=idx, columns=idx)
    .fillna(0, downcast='infer')
    .pipe(lambda x : x+x.values.T))

# LWED
relatedness_file_LWED = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/plink/king_related_LWED.kin0"
related_df_LWED = pd.read_csv(relatedness_file_LWED, sep = "\t")
idx = sorted(set(related_df_LWED["FID1"]).union(related_df_LWED["FID2"]))
matrix_LWED= (related_df_LWED.pivot(index="FID1", columns = "FID2", values="Kinship")
    .reindex(index=idx, columns=idx)
    .fillna(0, downcast='infer')
    .pipe(lambda x : x+x.values.T))




matrix_SIMP.to_csv("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/plink/relatedness_matrix_SIMP.tsv", sep="\t")
matrix_LWED.to_csv("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/plink/relatedness_matrix_LWED.tsv", sep="\t")


import seaborn as sns; sns.set()
import matplotlib.pyplot as plt

sns.heatmap(matrix)
ax = sns.clustermap(matrix_LWED, annot= True)


from matplotlib.colors import BoundaryNorm, ListedColormap, LinearSegmentedColormap
from matplotlib import colors

# attempt - works! kind of
#cmap = sns.cubehelix_palette(start=2.8, rot=.1, light = 0.9, n_colors=3)
#grid_kws = {'width_ratios': (0.9,0.1), 'wspace': .1}
#fig, (ax, cbar_ax) = plt.subplots(1,2, gridspec_kw=grid_kws)
#ax = sns.heatmap(matrix_SIMP, ax = ax, cbar_ax=cbar_ax, cmap=ListedColormap(cmap), \
#    linewidths=0.5, linecolor="lightgray", cbar_kws={'orientation':'vertical'})
#ax.set_xticklabels(ax.get_xticklabels(), rotation = 45)


# LWED
lwed_mask = np.zeros_like(matrix_LWED, dtype=np.bool)
lwed_mask[np.triu_indices_from(lwed_mask)] = True
lwed_mask[np.diag_indices_from(lwed_mask)] = False

cmap = sns.cubehelix_palette(start=2.8, rot=.1, light = 0.9, n_colors=3)
sns.set(rc = {"figure.figsize":(10,8)})
ax = sns.heatmap(matrix_LWED, cmap=cmap, mask=lwed_mask, annot= True, annot_kws={"size":9})
colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0,0.177,0.354])
colorbar.set_ticklabels([">= 2nd degree", "1st degree", "MZ twin/duplicate"])
ax.set_title("Pairwise kinship coefficients for L. weddelli")
ax.figure.tight_layout()
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/kinship_heatmap_LWED_aug.png", dpi=300)


# SIMP
simp_mask = np.zeros_like(matrix_SIMP, dtype=np.bool)
simp_mask[np.triu_indices_from(simp_mask)] = True
simp_mask[np.diag_indices_from(simp_mask)] = False

cmap = sns.cubehelix_palette(start=2.8, rot=.1, light = 0.9, n_colors=3)
sns.set(rc = {"figure.figsize":(10,8)})
ax = sns.heatmap(matrix_SIMP, cmap=cmap, mask=simp_mask, annot= True, annot_kws={"size":10})
colorbar = ax.collections[0].colorbar
colorbar.set_ticks([0,0.177,0.354])
colorbar.set_ticklabels([">= 2nd degree", "1st degree", "MZ twin/duplicate"])
ax.set_title("Pairwise kinship coefficients for S. imperator")
ax.figure.tight_layout()
plt.savefig("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/kinship_heatmap_SIMP.png", dpi=300)







matrix.to_csv("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/relatedness_matrix.tsv", sep="\t")

# okay so we know that animal 10 has samples (4,32,26) which all show up as duplicates in king.
# animal 73 has samples (29,1) (dups in KING) but seems to be highly related to sample 18 and 23 from individuals 466 and 1223 respectively. All four are from the same group! IC
# animal 182 has samples (34, 28, 6) (dups in KING) but seems to be highly related to sample 17 and 19 from inds 516 and 429 respectively. All five are from same group, Ji
# animal 203 has samples (30, 2) (dups in KING)
# animal 27 has samples (27, 33, 5) (dups in KING)
# animal 216 has samples (31, 3).

# we could potentially now compare repeated vs potential siblings. 
# let's first choose which of the repeated samples out of each rep set we will keep, ideally the one with the less xs in the sequence.
snp_seqs = "/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/dual_ind_test_1.tsv"
snp_seqs_df = pd.read_csv(snp_seqs, sep = "\t")
x_count= {}
for index,rows in snp_seqs_df.iterrows():
    ind_count = 0
    for pos in rows[1]:
        if pos =="N":
            ind_count = ind_count + 1
    x_count[rows[0]] = ind_count


# keep the follwing 10 individuals:
tamOpt-003_S3,tamOpt-006_S6,tamOpt-010_S9,tamOpt-012_S10,tamOpt-014_S12,tamOpt-015_S13,tamOpt-016_S14,tamOpt-020_S21,tamOpt-024_S29,tamOpt-035_S33

# bcftools view -Oz -s tamOpt-003_S3,tamOpt-006_S6,tamOpt-010_S9,tamOpt-012_S10,tamOpt-014_S12,tamOpt-015_S13,tamOpt-016_S14,tamOpt-020_S21,tamOpt-024_S29,tamOpt-035_S33 vcf_stage_2_snpnames.vcf.gz > vcf_stage_2_snpnames_rmrel.vcf.gz

# bcftools +fill-tags vcf_stage_2_snpnames_rmrel.vcf.gz > vcf_stage_2_snpnames_rmrel_recalculated.vcf
# bgzip -c vcf_stage_2_snpnames_rmrel_recalculated.vcf > vcf_stage_2_snpnames_rmrel_recalculated.vcf.gz
# tabix -p vcf vcf_stage_2_snpnames_rmrel_recalculated.vcf.gz

###########################
from pysam import VariantFile

vcf_in = VariantFile("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/vcf_stage_2_snpnames_rmrel_recalculated.vcf.gz")
vcf_in_LWED = VariantFile("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/vcf_stage_2_snpnames_LWED_unrel_recalculated.vcf.gz")

# pasting here since can't load pysam and matplotlib in the same session
dictionary = {'21': (0.5, 6, 2, 4), '22': (None, 0, 0, 0), '24': (0.5, 8, 6, 2), '25': (0.4375, 8, 4, 5), '29': (None, 0, 0, 0), '30': (None, 0, 0, 0), '31': (0.4375, 8, 4, 3), '32': (0.33333298563957214, 9, 2, 4), '34': (0.38888901472091675, 9, 2, 5), '38': (0.5, 9, 4, 5), '39': (0.38888901472091675, 9, 0, 7), '43': (None, 0, 0, 0), '44': (0.5, 3, 2, 1), '46': (None, 0, 0, 0), '47': (0.5, 8, 2, 6), '49': (0.4000000059604645, 5, 0, 4), '53': (0.38888901472091675, 9, 6, 5), '54': (0.4375, 8, 2, 5), '55': (0.2777779996395111, 9, 4, 1), '60': (0.41666701436042786, 6, 4, 3), '61': (None, 0, 0, 0), '62': (0.20000000298023224, 10, 0, 4), '64': (0.4285709857940674, 7, 4, 4), '65': (0.4444440007209778, 9, 4, 4), '69': (None, 0, 0, 0), '70': (0.0, 2, 0, 0), '74': (0.5, 9, 4, 5), '75': (0.5, 9, 6, 3), '76': (0.4444440007209778, 9, 4, 4), '78': (0.375, 8, 0, 6), '80': (0.4444440007209778, 9, 0, 8), '81': (0.38888901472091675, 9, 6, 5), '82': (0.30000001192092896, 5, 2, 1), '83': (0.4444440007209778, 9, 6, 4), '86': (0.5, 1, 0, 1), '88': (0.3125, 8, 2, 3), '89': (0.25, 4, 4, 2), '90': (0.30000001192092896, 10, 0, 6), '91': (0.5, 2, 2, 0), '92': (0.4444440007209778, 9, 4, 4), '94': (0.33333298563957214, 9, 4, 2), '98': (None, 0, 0, 0), '100': (0.2777779996395111, 9, 10, 3), '101': (0.5, 1, 0, 1), '103': (0.5, 2, 2, 0), '104': (0.38888901472091675, 9, 2, 5), '105': (0.0, 1, 2, 0), '107': (None, 0, 0, 0), '108': (0.5, 5, 2, 3), '109': (0.5, 6, 2, 4), '110': (0.5, 3, 2, 1), '111': (0.4285709857940674, 7, 6, 2), '112': (0.4000000059604645, 10, 4, 4), '113': (0.5, 8, 6, 2), '114': (0.375, 8, 2, 4), '117': (0.375, 4, 2, 3), '118': (None, 0, 0, 0), '119': (0.4375, 8, 6, 3), '120': (None, 0, 0, 0), '123': (0.4375, 8, 2, 5), '124': (0.4375, 8, 2, 5), '128': (0.5, 8, 6, 2), '129': (None, 0, 0, 0), '130': (None, 0, 0, 0), '131': (0.38888901472091675, 9, 6, 5), '132': (None, 0, 0, 0), '133': (0.375, 4, 2, 3), '134': (0.0, 2, 0, 0), '135': (0.1875, 8, 0, 3), '137': (None, 0, 0, 0), '139': (None, 0, 0, 0), '140': (0.4375, 8, 2, 7), '141': (0.33333298563957214, 9, 8, 4), '142': (0.0, 2, 0, 0), '144': (0.4375, 8, 2, 5), '145': (0.25, 6, 8, 1), '146': (None, 0, 0, 0), '147': (0.30000001192092896, 5, 2, 1), '148': (0.5, 1, 0, 1), '149': (0.38888901472091675, 9, 4, 3), '150': (0.125, 8, 0, 2), '152': (0.375, 8, 6, 4), '153': (0.25, 2, 0, 1), '154': (None, 0, 0, 0), '156': (0.33333298563957214, 3, 4, 0), '157': (None, 0, 0, 0), '158': (0.5, 5, 4, 1), '159': (None, 0, 0, 0), '160': (0.25, 8, 0, 4), '161': (0.3125, 8, 6, 5), '162': (0.4375, 8, 2, 5), '164': (0.25, 6, 2, 1), '165': (0.0, 2, 4, 0), '166': (None, 0, 0, 0), '167': (0.5, 4, 2, 2), '168': (0.5, 2, 0, 2), '169': (0.4285709857940674, 7, 2, 6), '171': (0.30000001192092896, 5, 2, 1), '173': (0.4285709857940674, 7, 4, 4), '176': (0.5, 7, 4, 3), '177': (0.375, 8, 0, 6), '178': (0.4444440007209778, 9, 4, 6), '179': (0.38888901472091675, 9, 2, 5), '181': (0.5, 9, 4, 5), '183': (0.4375, 8, 4, 5), '184': (0.33333298563957214, 6, 4, 4), '187': (0.41666701436042786, 6, 0, 5), '188': (None, 0, 0, 0), '189': (0.1875, 8, 0, 3), '190': (0.21428599953651428, 7, 8, 3), '192': (0.4285709857940674, 7, 6, 2), '193': (0.0, 2, 4, 0), '358': (0.4285709857940674, 7, 4, 4), '360': (0.33333298563957214, 9, 0, 6), '363': (0.4444440007209778, 9, 10, 0), '364': (0.4375, 8, 4, 3), '367': (0.2857140004634857, 7, 2, 2), '369': (0.4444440007209778, 9, 2, 6), '372': (0.5, 7, 2, 5), '373': (0.125, 8, 0, 2), '376': (0.125, 8, 0, 2), '378': (0.5, 8, 8, 0), '379': (0.4444440007209778, 9, 2, 6), '380': (0.07142859697341919, 7, 0, 1), '381': (0.41666701436042786, 6, 2, 3), '382': (0.38888901472091675, 9, 6, 5), '383': (None, 0, 0, 0), '384': (0.4375, 8, 4, 5), '385': (0.38888901472091675, 9, 4, 3), '388': (0.4444440007209778, 9, 6, 4), '389': (0.33333298563957214, 9, 2, 4), '392': (0.4444440007209778, 9, 2, 6), '393': (0.4444440007209778, 9, 8, 2), '394': (0.4375, 8, 6, 3), '395': (0.4444440007209778, 9, 8, 0), '396': (0.4375, 8, 4, 3), '399': (0.35714301466941833, 7, 0, 5), '400': (0.33333298563957214, 9, 6, 6), '401': (0.5, 8, 4, 4), '293': (0.0, 9, 18, 0), '296': (0.4444440007209778, 9, 10, 0), '297': (0.4444440007209778, 9, 10, 0), '299': (0.11111100018024445, 9, 0, 2), '302': (0.5, 9, 8, 1), '304': (0.4000000059604645, 5, 4, 2), '306': (0.33333298563957214, 9, 8, 4), '307': (0.20000000298023224, 5, 0, 2), '308': (0.4444440007209778, 9, 10, 0), '311': (0.4444440007209778, 9, 8, 0), '312': (0.0, 9, 18, 0), '313': (0.2777779996395111, 9, 0, 5), '315': (0.0, 9, 0, 0), '317': (0.16666699945926666, 9, 2, 1), '320': (0.25, 8, 0, 4), '321': (0.2777779996395111, 9, 0, 5), '322': (0.25, 8, 2, 2), '323': (0.4444440007209778, 9, 10, 0), '324': (0.2222220003604889, 9, 10, 4), '325': (0.20000000298023224, 5, 2, 0), '327': (0.25, 8, 2, 2), '328': (0.11111100018024445, 9, 0, 2), '329': (0.4444440007209778, 9, 10, 0), '331': (0.0, 9, 0, 0), '332': (0.30000001192092896, 5, 0, 3), '333': (0.4285709857940674, 7, 4, 4), '334': (0.05555560067296028, 9, 0, 1), '335': (0.4444440007209778, 9, 8, 0), '336': (0.2222220003604889, 9, 2, 2), '337': (0.2222220003604889, 9, 0, 4), '341': (0.4444440007209778, 9, 8, 0), '345': (0.0, 5, 0, 0), '346': (0.5, 8, 8, 0), '348': (0.20000000298023224, 5, 2, 0), '352': (0.4285709857940674, 7, 8, 0), '355': (0.11111100018024445, 9, 0, 2), '356': (0.4444440007209778, 9, 8, 0), '357': (0.38888901472091675, 9, 6, 1), '226': (0.4444440007209778, 9, 10, 0), '232': (0.4285709857940674, 7, 4, 2), '233': (0.25, 4, 0, 2), '235': (0.11111100018024445, 9, 0, 2), '237': (0.375, 4, 2, 3), '239': (0.33333298563957214, 9, 4, 2), '240': (0.4444440007209778, 9, 8, 0), '244': (0.4000000059604645, 5, 6, 0), '246': (0.4444440007209778, 9, 8, 0), '248': (0.0, 1, 0, 0), '249': (0.07142859697341919, 7, 0, 1), '250': (0.33333298563957214, 3, 2, 0), '254': (0.25, 4, 4, 2), '255': (0.5, 8, 8, 0), '257': (0.2777779996395111, 9, 4, 1), '259': (0.2777779996395111, 9, 2, 3), '263': (0.5, 9, 8, 1), '264': (0.0625, 8, 0, 1), '266': (0.375, 8, 6, 0), '268': (0.33333298563957214, 9, 4, 2), '270': (None, 0, 0, 0), '271': (0.05555560067296028, 9, 0, 1), '272': (0.4444440007209778, 9, 8, 0), '274': (0.125, 4, 0, 1), '275': (0.5, 8, 8, 0), '276': (0.0, 1, 2, 0), '278': (0.33333298563957214, 9, 4, 2), '279': (0.4444440007209778, 9, 10, 0), '281': (0.2777779996395111, 9, 4, 1), '282': (0.33333298563957214, 6, 2, 2), '283': (0.4444440007209778, 9, 10, 0), '284': (0.2857140004634857, 7, 2, 2), '285': (0.30000001192092896, 5, 2, 1), '286': (0.14285700023174286, 7, 0, 2), '288': (0.05555560067296028, 9, 0, 1), '289': (0.05555560067296028, 9, 0, 1), '290': (0.16666699945926666, 9, 2, 1)}
dictionary = {}
dictionary_LWED = {'21': (0.33333298563957214, 3, 0, 2), '22': (None, 0, 0, 0), '24': (0.375, 4, 2, 1), '25': (0.5, 4, 2, 2), '29': (None, 0, 0, 0), '30': (None, 0, 0, 0), '31': (0.375, 4, 2, 3), '32': (0.375, 4, 2, 1), '34': (0.375, 4, 0, 3), '38': (0.375, 4, 0, 3), '39': (0.375, 4, 0, 3), '43': (None, 0, 0, 0), '44': (0.5, 2, 2, 0), '46': (None, 0, 0, 0), '47': (0.375, 4, 2, 3), '49': (0.25, 2, 0, 1), '53': (0.5, 4, 2, 2), '54': (0.375, 4, 2, 3), '55': (0.375, 4, 4, 1), '60': (0.375, 4, 2, 3), '61': (None, 0, 0, 0), '62': (0.30000001192092896, 5, 0, 3), '64': (0.33333298563957214, 3, 2, 2), '65': (0.375, 4, 2, 3), '69': (None, 0, 0, 0), '70': (0.0, 1, 0, 0), '74': (0.5, 5, 2, 3), '75': (0.125, 4, 6, 1), '76': (0.5, 4, 2, 2), '78': (0.375, 4, 0, 3), '80': (0.375, 4, 0, 3), '81': (0.5, 4, 0, 4), '82': (0.0, 2, 0, 0), '83': (0.25, 4, 4, 2), '86': (0.5, 1, 0, 1), '88': (0.5, 3, 2, 1), '89': (0.25, 2, 2, 1), '90': (0.30000001192092896, 5, 0, 3), '91': (0.0, 1, 0, 0), '92': (0.4000000059604645, 5, 2, 2), '94': (0.375, 4, 2, 1), '98': (None, 0, 0, 0), '100': (0.5, 4, 2, 2), '101': (0.5, 1, 0, 1), '103': (0.0, 1, 2, 0), '104': (0.30000001192092896, 5, 0, 3), '105': (0.0, 1, 2, 0), '107': (None, 0, 0, 0), '108': (0.5, 4, 2, 2), '109': (0.33333298563957214, 3, 2, 2), '110': (0.5, 2, 2, 0), '111': (0.375, 4, 4, 1), '112': (0.20000000298023224, 5, 0, 2), '113': (0.16666699945926666, 3, 0, 1), '114': (0.375, 4, 2, 1), '117': (0.5, 2, 0, 2), '118': (None, 0, 0, 0), '119': (0.25, 4, 4, 2), '120': (None, 0, 0, 0), '123': (0.5, 4, 2, 2), '124': (0.25, 4, 0, 2), '128': (0.33333298563957214, 3, 2, 0), '129': (None, 0, 0, 0), '130': (None, 0, 0, 0), '131': (0.375, 4, 0, 3), '132': (None, 0, 0, 0), '133': (0.5, 2, 0, 2), '134': (0.0, 1, 0, 0), '135': (0.25, 4, 0, 2), '137': (None, 0, 0, 0), '139': (None, 0, 0, 0), '140': (0.5, 3, 0, 3), '141': (0.375, 4, 2, 3), '142': (None, 0, 0, 0), '144': (0.25, 4, 0, 2), '145': (0.375, 4, 4, 1), '146': (None, 0, 0, 0), '147': (0.0, 3, 0, 0), '148': (None, 0, 0, 0), '149': (0.25, 4, 4, 2), '150': (0.25, 4, 0, 2), '152': (0.375, 4, 4, 1), '153': (0.5, 1, 0, 1), '154': (None, 0, 0, 0), '156': (0.5, 2, 2, 0), '157': (None, 0, 0, 0), '158': (0.16666699945926666, 3, 4, 1), '159': (None, 0, 0, 0), '160': (0.125, 4, 0, 1), '161': (0.0, 3, 6, 0), '162': (0.25, 4, 0, 2), '164': (0.33333298563957214, 3, 2, 0), '165': (0.0, 2, 4, 0), '166': (None, 0, 0, 0), '167': (0.33333298563957214, 3, 2, 2), '168': (None, 0, 0, 0), '169': (0.5, 2, 0, 2), '171': (0.25, 4, 2, 0), '173': (0.5, 4, 2, 2), '176': (0.33333298563957214, 3, 0, 2), '177': (0.25, 4, 0, 2), '178': (0.375, 4, 0, 3), '179': (0.25, 4, 0, 2), '181': (0.25, 4, 0, 2), '183': (0.5, 4, 0, 4), '184': (0.16666699945926666, 3, 4, 1), '187': (0.33333298563957214, 3, 0, 2), '188': (None, 0, 0, 0), '189': (0.25, 4, 0, 2), '190': (0.16666699945926666, 3, 4, 1), '192': (0.375, 4, 2, 1), '193': (0.0, 1, 2, 0), '358': (0.5, 3, 0, 3), '360': (0.5, 4, 0, 4), '363': (0.0, 4, 0, 0), '364': (0.25, 4, 0, 2), '367': (0.16666699945926666, 3, 0, 1), '369': (0.25, 4, 0, 2), '372': (0.25, 2, 0, 1), '373': (0.16666699945926666, 3, 0, 1), '376': (0.0, 3, 0, 0), '378': (0.0, 4, 0, 0), '379': (0.375, 4, 2, 3), '380': (0.25, 2, 0, 1), '381': (0.33333298563957214, 3, 0, 2), '382': (0.125, 4, 6, 1), '383': (None, 0, 0, 0), '384': (0.16666699945926666, 3, 4, 1), '385': (0.0, 4, 0, 0), '388': (0.375, 4, 0, 3), '389': (0.0, 4, 0, 0), '392': (0.5, 4, 2, 2), '393': (0.0, 4, 8, 0), '394': (0.16666699945926666, 3, 4, 1), '395': (0.0, 4, 8, 0), '396': (0.25, 4, 4, 2), '399': (0.33333298563957214, 3, 0, 2), '400': (0.25, 4, 4, 2), '401': (0.33333298563957214, 3, 2, 0), '293': (0.0, 4, 8, 0), '296': (0.0, 4, 0, 0), '297': (0.0, 4, 0, 0), '299': (0.0, 4, 0, 0), '302': (0.0, 4, 0, 0), '304': (None, 0, 0, 0), '306': (0.375, 4, 0, 3), '307': (None, 0, 0, 0), '308': (0.0, 4, 0, 0), '311': (0.0, 4, 8, 0), '312': (0.0, 4, 8, 0), '313': (0.0, 4, 0, 0), '315': (0.0, 4, 0, 0), '317': (0.0, 4, 0, 0), '320': (0.0, 4, 0, 0), '321': (0.125, 4, 0, 1), '322': (0.5, 4, 2, 2), '323': (0.0, 4, 0, 0), '324': (0.0, 4, 8, 0), '325': (None, 0, 0, 0), '327': (0.0, 4, 0, 0), '328': (0.0, 4, 0, 0), '329': (0.0, 4, 0, 0), '331': (0.0, 4, 0, 0), '332': (None, 0, 0, 0), '333': (0.0, 2, 4, 0), '334': (0.0, 4, 0, 0), '335': (0.0, 4, 8, 0), '336': (0.0, 4, 0, 0), '337': (0.0, 4, 0, 0), '341': (0.0, 4, 8, 0), '345': (None, 0, 0, 0), '346': (0.0, 4, 0, 0), '348': (0.0, 3, 0, 0), '352': (0.0, 3, 0, 0), '355': (0.0, 4, 0, 0), '356': (0.0, 4, 8, 0), '357': (0.0, 4, 0, 0), '226': (0.0, 4, 0, 0), '232': (0.25, 4, 4, 2), '233': (0.33333298563957214, 3, 0, 2), '235': (0.25, 4, 0, 2), '237': (0.375, 4, 2, 3), '239': (0.25, 4, 4, 2), '240': (0.0, 4, 8, 0), '244': (0.0, 3, 6, 0), '246': (0.0, 4, 8, 0), '248': (0.0, 1, 0, 0), '249': (0.16666699945926666, 3, 0, 1), '250': (0.0, 2, 0, 0), '254': (0.25, 4, 4, 2), '255': (0.0, 4, 0, 0), '257': (0.375, 4, 4, 1), '259': (0.375, 4, 2, 3), '263': (0.0, 4, 8, 0), '264': (0.16666699945926666, 3, 0, 1), '266': (0.0, 3, 6, 0), '268': (0.25, 4, 4, 2), '270': (None, 0, 0, 0), '271': (0.10000000149011612, 5, 0, 1), '272': (0.0, 4, 8, 0), '274': (0.125, 4, 0, 1), '275': (0.0, 4, 0, 0), '276': (0.0, 1, 2, 0), '278': (0.25, 4, 4, 2), '279': (0.0, 4, 0, 0), '281': (0.375, 4, 4, 1), '282': (0.33333298563957214, 3, 2, 2), '283': (0.0, 4, 0, 0), '284': (0.5, 4, 2, 2), '285': (0.5, 3, 2, 1), '286': (0.25, 4, 0, 2), '288': (0.125, 4, 0, 1), '289': (0.125, 4, 0, 1), '290': (0.375, 4, 2, 1)}
dictionary_LWED = {}

for rec in vcf_in_LWED.fetch():
    print(rec.info.keys())
    snp= rec.id
    maf = rec.info["MAF"]
    ns = rec.info["NS"]
    ac_hom = rec.info["AC_Hom"][0]
    ac_het = rec.info["AC_Het"][0]
    dictionary_LWED[snp] = maf , ns, ac_hom , ac_het

# let's look at the maf distribution
list_mafs = []
for snp in dictionary_LWED:
    if dictionary_LWED[snp][0] is not None:
        list_mafs.append(dictionary_LWED[snp][0])

import seaborn as sns; sns.set()
import matplotlib.pyplot as plt

sns.boxplot(list_mafs)

# let's have a closer look at the sites that have a low MAF, let's say, below or equal to 0.29
low_maf= []
for snp in dictionary_LWED:
    if dictionary_LWED[snp][0] is not None and dictionary_LWED[snp][0] != 0:
        if float(dictionary_LWED[snp][0]) <= 0.29:
            low_maf.append(snp)

# low_maf = >>> low_maf
# ['55', '62', '89', '100', '135', '145', '150', '153', '160', '164', '189', '190', '367', '373', '376', '380', '299', '307', '313', '317', '320', '321', '322', '324', '325', '327', '328', '334', '336', '337', '348', '355', '233', '235', '249', '254', '257', '259', '264', '271', '274', '281', '284', '286', '288', '289', '290']
# # len(low_maf) is 47

no_maf = []
maf_obtained = []
for snp in dictionary_LWED:
    if dictionary_LWED[snp][0] is None:
        no_maf.append(snp)
    elif dictionary_LWED[snp][0] is not None:
        maf_obtained.append(snp)

#len(no_maf) is 32
#['22', '26', '29', '30', '37', '42', '43', '46', '61', '67', '69', '98', '107', '118', '120', '129', '130', '132', '137', '138', '139', '146', '154', '157', '159', '166', '174', '188', '361', '383', '326', '270']

# len(maf_obtained) is 190
# ['21', '24', '25', '31', '32', '34', '38', '39', '44', '47', '49', '53', '54', '55', '60', '62', '64', '65', '70', '74', '75', '76', '78', '80', '81', '82', '83', '86', '88', '89', '90', '91', '92', '94', '100', '101', '103', '104', '105', '108', '109', '110', '111', '112', '113', '114', '117', '119', '123', '124', '128', '131', '133', '134', '135', '140', '141', '142', '144', '145', '147', '148', '149', '150', '152', '153', '156', '158', '160', '161', '162', '164', '165', '167', '168', '169', '171', '173', '176', '177', '178', '179', '181', '183', '184', '187', '189', '190', '192', '193', '358', '360', '363', '364', '367', '369', '372', '373', '376', '378', '379', '380', '381', '382', '384', '385', '388', '389', '392', '393', '394', '395', '396', '399', '400', '401', '293', '296', '297', '299', '302', '304', '306', '307', '308', '311', '312', '313', '315', '317', '320', '321', '322', '323', '324', '325', '327', '328', '329', '331', '332', '333', '334', '335', '336', '337', '341', '345', '346', '348', '352', '355', '356', '357', '226', '232', '233', '235', '237', '239', '240', '244', '246', '248', '249', '250', '254', '255', '257', '259', '263', '264', '266', '268', '271', '272', '274', '275', '276', '278', '279', '281', '282', '283', '284', '285', '286', '288', '289', '290']

# another way to parse the vcf since I want the actual genotypes, maf shouldn't be our only indicator!
############################

from cyvcf2 import VCF # if it doesnt work go to glob and run from there

variant_dict = {}
for variant in VCF("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/vcf_stage_2_snpnames_rmrel_recalculated.vcf.gz"):
    hom_ref = 0
    hom_alt = 0
    het = 0
    no_call= 0
    ns = 0  # number of samples with data
    for i in range(9,19):
        if str(variant).split()[i] == './.':
            no_call = no_call + 1
        elif str(variant).split()[i] == '0/0':
            hom_ref = hom_ref + 1
            ns = ns + 1
        elif (str(variant).split()[i] == '1/1') or (str(variant).split()[i] == "2/2"):
            hom_alt = hom_alt + 1
            ns = ns + 1
        elif (str(variant).split()[i] == '0/1') or (str(variant).split()[i] == "1/2"):
            het = het + 1
            ns = ns + 1
        else:
            print(str(variant).split())
    variant_dict[variant.ID] = hom_ref, hom_alt, het, no_call, ns

# LWED and SIMP
variant_dict = {}
for variant in VCF("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/vcf_stage_2_snpnames_LWED_unrel_recalculated.vcf.gz"):
    hom_ref = 0
    hom_alt = 0
    het = 0
    no_call= 0
    ns = 0  # number of samples with data
    for i in range(9,13):
        if str(variant).split()[i] == './.':
            no_call = no_call + 1
        elif str(variant).split()[i] == '0/0':
            hom_ref = hom_ref + 1
            ns = ns + 1
        elif (str(variant).split()[i] == '1/1') or (str(variant).split()[i] == "2/2"):
            hom_alt = hom_alt + 1
            ns = ns + 1
        elif (str(variant).split()[i] == '0/1') or (str(variant).split()[i] == "1/2"):
            het = het + 1
            ns = ns + 1
        else:
            print(str(variant).split())
    variant_dict[variant.ID] = hom_ref, hom_alt, het, no_call, ns




# looking at what the variants with low MAF look like:

low_maf = ['55', '62', '89', '100', '135', '145', '150', '153', '160', '164', '189', '190', '367', '373', '376', '380', '299', '307', '313', '317', '320', '321', '322', '324', '325', '327', '328', '334', '336', '337', '348', '355', '233', '235', '249', '254', '257', '259', '264', '271', '274', '281', '284', '286', '288', '289', '290']
low_maf_dict = {}
for variant in variant_dict:
    if variant in low_maf:
        low_maf_dict[variant] = variant_dict[variant]

# want to make a table with this so here we go:
for variant in low_maf_dict:
    print(str(variant) + "," + str(low_maf_dict[variant][0]) + "," + str(low_maf_dict[variant][1]) + "," + str(low_maf_dict[variant][2]) + "," + str(low_maf_dict[variant][3]) + "," + str(low_maf_dict[variant][4]))



# ok, calculated! 
# again, let's remove all sites for which no calls equal 0
no_calls = []
for snp in variant_dict:
    if variant_dict[snp][4] == 0:
        no_calls.append(snp)

>>> len(no_calls)
24

more_than_point3_calls = []
less_than_3 = []
more_than_3 = {}
for snp in variant_dict:
    if variant not in no_calls:
        if variant_dict[snp][4] >= 3:
            more_than_point3_calls.append(snp)
            more_than_3[snp]=variant_dict[snp]
        else:
            less_than_3.append(snp)
# we have 175 variants producing genotypes in more than 30% of the unrelated samples!

# let's see how many variants are presenting ALL of their genotypes as hets:
for snp in variant_dict:
    if variant_dict[snp][4] > 0 and variant_dict[snp][4] == variant_dict[snp][2]:
        print(snp, variant_dict[snp])

# so we have these sites:
86 (0, 0, 1, 9, 1)
101 (0, 0, 1, 9, 1)
148 (0, 0, 1, 9, 1)
168 (0, 0, 2, 8, 2)

# I don't see a need to remove these other than the low amount of genotypes produced.

# let's see how many variants are presenting ALL their genotypes as homs for ref:
for snp in variant_dict:
    if variant_dict[snp][4] > 0 and variant_dict[snp][4] == variant_dict[snp][0]:
        print(snp, variant_dict[snp])

# we got these
70 (2, 0, 0, 8, 2)
134 (2, 0, 0, 8, 2)
142 (2, 0, 0, 8, 2)
315 (9, 0, 0, 1, 9)*
331 (9, 0, 0, 1, 9)*
345 (5, 0, 0, 5, 5)*
248 (1, 0, 0, 9, 1)

# again, low amount of genotypes called but some of them a bit concerning.

# let's see how many variants are presenting ALL their genotypes as homs for alt:
for snp in variant_dict:
    if variant_dict[snp][4] > 0 and variant_dict[snp][4] == variant_dict[snp][1]:
        print(snp, variant_dict[snp])

# we got these two:
105 (0, 1, 0, 9, 1)
165 (0, 2, 0, 8, 2)
193 (0, 2, 0, 8, 2)
293 (0, 9, 0, 1, 9)*
312 (0, 9, 0, 1, 9)*
276 (0, 1, 0, 9, 1)

# again, low amount of reads but some a bit concerning.

# let's see how many of our variants are presenting the same genotype 90% of the time: .
# HOM REF
for variant in variant_dict: 
    if variant_dict[variant][4] > 0 and variant_dict[variant][0] >= variant_dict[variant][4]*.8: 
        print(variant, variant_dict[variant])

#90% 
SAME AS 100%

#80%:
380 (6, 0, 1, 3, 7)*
325 (4, 1, 0, 5, 5)
334 (8, 0, 1, 1, 9)*
249 (6, 0, 1, 3, 7)*
264 (7, 0, 1, 2, 8)*
271 (8, 0, 1, 1, 9)*
288 (8, 0, 1, 1, 9)*
289 (8, 0, 1, 1, 9)*


# HOM ALT
for variant in variant_dict:
    if variant_dict[variant][4] > 0 and variant_dict[variant][1] >= variant_dict[variant][4]*.8: # 90% nothing
        print(variant, variant_dict[variant])
#90%
SAME AS 100%

#80%
SAME AS 100%


# HET
for variant in variant_dict:
    if variant_dict[variant][4] > 0 and variant_dict[variant][2] >= variant_dict[variant][4]*.8: 
        print(variant, variant_dict[variant])

#90%
SAME AS 100%

#80%
49 (1, 0, 4, 5, 5)
80 (1, 0, 8, 1, 9)*
140 (0, 1, 7, 2, 8)*
169 (0, 1, 6, 3, 7)*
187 (1, 0, 5, 4, 6)*


# let's put some of the worrysome variants in a list:
worrysome = ["315", "331", "345", "293", "312", "380", "334", "249", "264", "271", "288", "289", "80", "140", "169", "187"]
worrysome = ["315", "331", "345", "293", "312", "380", "334", "249", "264", "271", "288", "289", "80", "140", "169", "187"]

# and we know the list of variants that we added earlier. 

rescued = ['393', '226', '240', '244', '246', '248', '250', '255', '263', '266', '272', '275', '279', '283', '293', '296', '297', '307', '308', '311', '312', '315', '316', '317', '318', '322', '323', '325', '329', '331', '334', '335', '341', '345', '352', '356']

for i in rescued:
    if i in worrysome:
        print(i)



# SAME THING BUT WANTED TO USE ALL INDIVIDUALS FOR THE NO CALLS FILTER!!!

#bcftools +fill-tags vcf_stage_2_snpnames.vcf.gz > vcf_stage_2_snpnames_recalculated.vcf
#bgzip -c vcf_stage_2_snpnames_recalculated.vcf > vcf_stage_2_snpnames_recalculated.vcf.gz
#tabix -p vcf vcf_stage_2_snpnames_recalculated.vcf.gz 

variant_dict = {}
for variant in VCF("/proj/proj_name/nobackup/SAM/genotyping/fourth_geno_run/vcf_stage_2_snpnames_recalculated.vcf.gz"):
    hom_ref = 0
    hom_alt = 0
    het = 0
    no_call= 0
    ac = 0
    for i in range(9,42):
        if str(variant).split()[i] == './.':
            no_call = no_call + 1
        elif str(variant).split()[i] == '0/0':
            hom_ref = hom_ref + 1
            ac = ac + 1
        elif str(variant).split()[i] == '1/1' or str(variant).split()[i] == "2/2":
            hom_alt = hom_alt + 1
            ac = ac + 1
        elif str(variant).split()[i] == '0/1' or str(variant).split()[i] == "1/2":
            het = het + 1
            ac = ac + 1
        else:
            print(str(variant).split())
    variant_dict[variant.ID] = hom_ref, hom_alt, het, no_call, ac

# ok, calculated! 
# again, let's remove all sites for which no calls equal 0
no_calls_2 = []
for snp in variant_dict:
    if variant_dict[snp][4] == 0:
        no_calls_2.append(snp)

>>> len(no_calls)
8

# now len(no_calls) is 2

# THIS IS GOOD, WE SHOULDN'T REMOVE THOSE 28, JUST THESE 9



>>>
# great, this is the same thing we had obtained earlier.

more_than_point3_calls = []
less_than_3 = []
more_than_3 = {}
zero = []
for snp in variant_dict:
    if variant_dict[snp][4] == 0:
        zero.append(snp)
    elif variant_dict[snp][4] >= 3:
        more_than_point3_calls.append(snp)
        more_than_3[snp]=variant_dict[snp]
    else:
        less_than_3.append(snp)
