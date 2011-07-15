#include "DataFormats/Common/interface/Wrapper.h"

//Add includes for your classes here
#include "MyDataFormats/RA2bDataFormats/interface/SimpleEt.h"

namespace {
   struct MyDataFormats_RA2bDataFormats {
     SimpleEt etdummy0;
     edm::Wrapper<SimpleEt> etdummy1;
   };
}
