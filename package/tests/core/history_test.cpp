#include "libcoherent4cast/history.hpp"
#include <cassert>
int main(){
 auto x=coherent4cast::compile_history({1,1,1,2},{0,1,3,0},{1,2,4,1},{1,2,0,0},2,2,.5);
 assert(x.counts[0]==0&&x.counts[2]==1&&x.counts[4]==1&&x.counts[5]==1&&x.counts[6]==0);
 assert(std::abs(x.recency[4]-std::exp(-1.0))<1e-14);
 assert(std::abs(x.recency[5]-std::exp(-.5))<1e-14);
 assert(x.gap[4]==2&&x.gap[5]==1&&std::isnan(x.gap[6]));
 bool failed=false;try{coherent4cast::compile_history({1,1},{0,.5},{1,2},{0,0},1,2,0);}catch(const std::invalid_argument&){failed=true;}assert(failed);
}
