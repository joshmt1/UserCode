#ifndef MyDataFormats_RA2bDataFormats_SimpleEt_h
#define MyDataFormats_RA2bDataFormats_SimpleEt_h


// a simple class
struct SimpleEt
{
  explicit SimpleEt(float v):pt(v),phi(v) { }
  SimpleEt():pt(-99),phi(-99) { }
  float pt,phi;
};


//typedef std::vector<SimpleEt> SampleCollection;

#endif
