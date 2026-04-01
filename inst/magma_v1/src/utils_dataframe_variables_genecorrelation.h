/** Copyright (C) 2014-2020 by Christiaan de Leeuw (CTG Lab, Vrije Universiteit Amsterdam), All Rights Reserved **/

#ifndef UTILS_DATAFRAME_VARIABLES_GENECORRELATION_H
#define UTILS_DATAFRAME_VARIABLES_GENECORRELATION_H

class GeneCorrelationVariable : public JaggedMultiVariable<PosCorrelationVariable<double> > {
  void reorder_insert(BaseBuffer<BaseBuffer<double>*>& corrs, int row, int col, double value);
  void reorder_pad(BaseBuffer<double>*& buff, int required);
public:
  GeneCorrelationVariable() : JaggedMultiVariable<PosCorrelationVariable<double> >(false) {}

  ///assumes no duplications in vector use
  void shrink(vector<int>& use);
  void reorder(vector<int>& use);
  
  void set_length(int iid, int length);
  
  JaggedVariableOutput<PosCorrelationVariable<double> >* get_output() {return !output_mask ? JaggedMultiVariable<PosCorrelationVariable<double> >::get_output() : 0;}
};


#endif /** UTILS_DATAFRAME_VARIABLES_GENECORRELATION_H */
