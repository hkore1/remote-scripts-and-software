#!/bin/bash

# From https://github.com/mossmatters/phyloscripts/tree/master/alleles_workflow

### Unhash these if running without a wrapper ###
#module load samtools/1.16.1
#module load bcftools/1.15
#module load parallel/20210322
#set -eo pipefail


#Script to prepare phased haplotype sequences for each for one sample. 

prefix=$1
genelist=$2

cd $prefix

# Create genelist.txt for target genes (currently configured for Angiosperm353 genes)
if [[ "$genelist" == "Angiosperms353" ]];
then
  rm -f $genelist
  touch $genelist
  echo -e "4471\n4527\n4691\n4724\n4744\n4757\n4793\n4796\n4802\n4806\n4848\n4889\n4890\n4893\n4932\n4942\n4951\n4954\n4989\n4992\n5018\n5032\n5034\n5038\n5064\n5090\n5104\n5116\n5123\n5131\n5138\n5162\n5163\n5168\n5177\n5188\n5200\n5206\n5220\n5257\n5260\n5264\n5271\n5273\n5280\n5296\n5299\n5304\n5318\n5326\n5328\n5333\n5335\n5339\n5343\n5347\n5348\n5354\n5355\n5357\n5366\n5398\n5404\n5406\n5421\n5422\n5426\n5427\n5428\n5430\n5434\n5449\n5454\n5460\n5463\n5464\n5469\n5477\n5489\n5502\n5513\n5528\n5531\n5536\n5551\n5554\n5562\n5578\n5594\n5596\n5599\n5614\n5620\n5634\n5639\n5642\n5644\n5656\n5660\n5664\n5670\n5699\n5702\n5703\n5716\n5721\n5733\n5744\n5770\n5772\n5791\n5802\n5815\n5816\n5821\n5822\n5840\n5841\n5842\n5843\n5849\n5853\n5857\n5858\n5859\n5865\n5866\n5870\n5893\n5894\n5899\n5910\n5913\n5918\n5919\n5921\n5922\n5926\n5933\n5936\n5940\n5941\n5942\n5943\n5944\n5945\n5949\n5950\n5958\n5960\n5968\n5974\n5977\n5980\n5981\n5990\n6000\n6003\n6004\n6016\n6026\n6029\n6034\n6036\n6038\n6041\n6048\n6050\n6051\n6056\n6064\n6068\n6072\n6098\n6110\n6114\n6119\n6128\n6130\n6139\n6148\n6150\n6164\n6175\n6176\n6198\n6216\n6221\n6226\n6227\n6238\n6258\n6265\n6270\n6274\n6282\n6284\n6295\n6298\n6299\n6303\n6318\n6320\n6363\n6366\n6373\n6376\n6378\n6379\n6383\n6384\n6387\n6389\n6393\n6398\n6401\n6404\n6405\n6406\n6407\n6412\n6420\n6430\n6432\n6439\n6447\n6448\n6449\n6450\n6454\n6457\n6458\n6459\n6460\n6462\n6483\n6487\n6488\n6492\n6494\n6496\n6498\n6500\n6506\n6507\n6514\n6526\n6527\n6528\n6531\n6532\n6533\n6538\n6540\n6544\n6550\n6552\n6557\n6559\n6563\n6565\n6570\n6572\n6601\n6620\n6631\n6636\n6639\n6641\n6649\n6652\n6660\n6667\n6679\n6685\n6689\n6705\n6713\n6717\n6732\n6733\n6738\n6746\n6779\n6780\n6782\n6785\n6791\n6792\n6797\n6825\n6848\n6854\n6859\n6860\n6864\n6865\n6875\n6882\n6883\n6886\n6893\n6909\n6913\n6914\n6924\n6933\n6946\n6947\n6954\n6955\n6958\n6961\n6962\n6968\n6969\n6977\n6978\n6979\n6992\n6995\n7013\n7021\n7024\n7028\n7029\n7067\n7111\n7128\n7135\n7136\n7141\n7174\n7194\n7241\n7273\n7279\n7296\n7313\n7324\n7325\n7331\n7333\n7336\n7361\n7363\n7367\n7371\n7572\n7577\n7583\n7602\n7628" > $genelist
fi

# Run bcftools to extract sequences

bgzip -c $prefix.supercontigs.fasta.snps.whatshap.vcf > $prefix.supercontigs.fasta.snps.whatshap.vcf.gz
tabix $prefix.supercontigs.fasta.snps.whatshap.vcf.gz
mkdir -p phased_bcftools
rm -f phased_bcftools/*

parallel "samtools faidx $prefix.supercontigs.fasta {1}---$prefix | bcftools consensus -H 1 $prefix.supercontigs.fasta.snps.whatshap.vcf.gz > phased_bcftools/$prefix-{1}_h1.phased.fasta" :::: $genelist 
parallel "samtools faidx $prefix.supercontigs.fasta {1}---$prefix | bcftools consensus -H 2 $prefix.supercontigs.fasta.snps.whatshap.vcf.gz >> phased_bcftools/$prefix-{1}_h2.phased.fasta" :::: $genelist 

wait
rm -f phased_bcftools/$prefix-.phased.fasta
cat phased_bcftools/*_h1.phased.fasta > phased_bcftools/"$prefix"_h1.fasta
cat phased_bcftools/*_h2.phased.fasta > phased_bcftools/"$prefix"_h2.fasta
rm -f phased_bcftools/*_h1.phased.fasta
rm -f phased_bcftools/*_h2.phased.fasta

cd ..
