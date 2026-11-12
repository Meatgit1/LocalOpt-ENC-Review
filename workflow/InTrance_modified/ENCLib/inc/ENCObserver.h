// Opt-in count-bank observer. No evaluator, controller, or random-generator calls.
#ifndef INTRANCE_ENC_OBSERVER_H
#define INTRANCE_ENC_OBSERVER_H

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <stdint.h>
#include <sys/stat.h>
#include <vector>
#ifdef _WIN32
#include <io.h>
#else
#include <dirent.h>
#include <unistd.h>
#endif

namespace ENCObservation {
inline void require(bool ok, const char* detail) {
  if (!ok) throw std::runtime_error(std::string("ENC_OBSERVER_REFUSE_") + detail);
}
inline std::string trim(const std::string& s) {
  const std::string::size_type first=s.find_first_not_of(" \t\r\n");
  if (first==std::string::npos) return "";
  return s.substr(first,s.find_last_not_of(" \t\r\n")-first+1);
}
inline unsigned long number(const std::string& text) {
  require(!text.empty() && text.find_first_not_of("0123456789")==std::string::npos,"INTEGER");
  errno=0; char* end=0; const unsigned long value=std::strtoul(text.c_str(),&end,10);
  require(errno!=ERANGE && end && *end==0 && value<=2147483647UL,"INTEGER_RANGE");
  return value;
}
inline bool finite(double x) {
  return x==x && x<=std::numeric_limits<double>::max() && x>=-std::numeric_limits<double>::max();
}
inline void emptyDirectory(const std::string& path) {
#ifdef _WIN32
  struct _stat st;
  require(_stat(path.c_str(),&st)==0 && (st.st_mode&_S_IFDIR)!=0,"OUTPUT_DIRECTORY");
  struct _finddata_t item;
  intptr_t handle=_findfirst((path+"/*").c_str(),&item);
  if (handle==-1) { require(errno==ENOENT,"OUTPUT_LIST"); return; }
  bool empty=true;
  do { if(std::string(item.name)!="." && std::string(item.name)!="..") empty=false; }
  while(_findnext(handle,&item)==0);
  _findclose(handle); require(empty,"OUTPUT_NOT_EMPTY");
#else
  struct stat st;
  require(lstat(path.c_str(),&st)==0 && S_ISDIR(st.st_mode),"OUTPUT_DIRECTORY");
  DIR* dir=opendir(path.c_str()); require(dir!=0,"OUTPUT_LIST");
  bool empty=true; struct dirent* item;
  while((item=readdir(dir))!=0) {
    if(std::string(item->d_name)!="." && std::string(item->d_name)!="..") empty=false;
  }
  closedir(dir); require(empty,"OUTPUT_NOT_EMPTY");
#endif
}
// Exclusive creation preserves previous attempts and incomplete evidence.
inline FILE* exclusive(const std::string& path) {
#ifdef _WIN32
  int fd=_open(path.c_str(),_O_CREAT|_O_EXCL|_O_WRONLY|_O_BINARY,_S_IREAD|_S_IWRITE);
  require(fd>=0,"FILE_EXISTS_OR_CREATE"); FILE* file=_fdopen(fd,"wb");
  if (!file) _close(fd);
#else
  int fd=open(path.c_str(),O_CREAT|O_EXCL|O_WRONLY,0600);
  require(fd>=0,"FILE_EXISTS_OR_CREATE"); FILE* file=fdopen(fd,"wb");
  if (!file) close(fd);
#endif
  require(file!=0,"FILE_OPEN"); return file;
}
inline void append(FILE* file,const std::string& text) {
  require(file && fwrite(text.data(),1,text.size(),file)==text.size() && fflush(file)==0,"FILE_WRITE");
}
inline void write(const std::string& path,const std::string& text) {
  FILE* f=exclusive(path);
  try { append(f,text); } catch(...) { fclose(f); throw; }
  require(fclose(f)==0,"FILE_CLOSE");
}
inline std::string controller(const std::vector<double>& values) {
  require(values.size()==881,"PARAMETER_COUNT");
  std::ostringstream s; s.imbue(std::locale::classic()); s<<std::setprecision(17);
  s<<values.size()<<"\n";
  for(size_t i=0;i<values.size();++i) { require(finite(values[i]),"NONFINITE_PARAMETER"); s<<values[i]<<"\n"; }
  return s.str();
}
inline std::string digest(const std::string& bytes) {
  unsigned long h=2166136261UL;
  for(size_t i=0;i<bytes.size();++i) h=((h^static_cast<unsigned char>(bytes[i]))*16777619UL)&0xffffffffUL;
  std::ostringstream s; s<<std::hex<<std::setw(8)<<std::setfill('0')<<h; return s.str();
}

struct Config {
  unsigned long seed,limit;
  bool capture;
  std::vector<unsigned long> checkpoints;
  std::string output;
  static Config load(const std::string& path) {
    require(!path.empty() && path.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.-")==std::string::npos && path!="." && path!="..","CONFIG_LOCAL_FILENAME");
    std::ifstream f(path.c_str()); require(f.is_open(),"CONFIG_OPEN");
    std::map<std::string,std::string> fields; std::string line;
    while(std::getline(f,line)) {
      line=trim(line); if(line.empty()) continue;
      size_t eq=line.find('='); require(eq!=std::string::npos && line.find('=',eq+1)==std::string::npos,"CONFIG_LINE");
      std::string key=trim(line.substr(0,eq)),value=trim(line.substr(eq+1));
      require(fields.insert(std::make_pair(key,value)).second,"DUPLICATE_KEY");
    }
    require(f.eof(),"CONFIG_READ");
    const char* keys[]={"schema","seed","eval_limit","capture","checkpoints","output_directory"};
    require(fields.size()==6,"CONFIG_KEYS");
    for(size_t i=0;i<6;++i) require(fields.find(keys[i])!=fields.end(),"CONFIG_KEYS");
    require(fields["schema"]=="1","SCHEMA");
    Config c; c.seed=number(fields["seed"]); c.limit=number(fields["eval_limit"]);
    require(c.limit>=40,"LIMIT_TOO_SMALL");
    require(fields["capture"]=="true" || fields["capture"]=="false","CAPTURE_BOOL");
    c.capture=fields["capture"]=="true"; c.output=fields["output_directory"];
    require(c.output=="enc_observer","OUTPUT_PATH");
    std::string marks=fields["checkpoints"]; require(!marks.empty() && marks[marks.size()-1]!=',',"CHECKPOINTS");
    std::istringstream list(marks); std::string mark;
    while(std::getline(list,mark,',')) {
      unsigned long n=number(trim(mark));
      require(n>0 && n<=c.limit && (c.checkpoints.empty() || n>c.checkpoints.back()),"CHECKPOINT_ORDER");
      c.checkpoints.push_back(n);
    }
    return c;
  }
};

class Observer {
  Config config_;
  unsigned long completed_,populations_,pairs_,reserved_,next_,last_boundary_;
  FILE* trace_; FILE* timing_;
  std::clock_t start_;
  std::string stop_;
  std::vector<unsigned long> actual_;
  bool finished_;
  Observer(const Observer&); Observer& operator=(const Observer&);
  std::string metadata(unsigned long requested,const char* reason,double fitness,
                       long scan,long epoch,double hc,double mutation,const std::string& hash) const {
    std::ostringstream s; s.imbue(std::locale::classic()); s<<std::setprecision(17);
    s<<"{\n\"schema\":1,\"seed\":"<<config_.seed<<",\"reason\":\""<<reason<<"\","
     <<"\"requested_evaluations\":"<<requested<<",\"completed_evaluations\":"<<completed_
     <<",\"completed_populations\":"<<populations_<<",\"completed_pairs\":"<<pairs_
     <<",\"scan_epoch\":"<<scan<<",\"refinement_epoch\":"<<epoch<<",\"hypercube\":"<<hc
     <<",\"mutation_scale\":"<<mutation<<",\"best_fitness\":"<<fitness
     <<",\"parameter_count\":881,\"controller_fnv1a32\":\""<<hash<<"\"}\n";
    return s.str();
  }
public:
  explicit Observer(const Config& config): config_(config),completed_(0),populations_(0),pairs_(0),reserved_(0),next_(0),last_boundary_(0),trace_(0),timing_(0),start_(std::clock()),finished_(false) {
    emptyDirectory(config_.output);
    trace_=exclusive(config_.output+"/trace.csv");
    try {
      timing_=exclusive(config_.output+"/timing.csv");
      append(trace_,"completed_evaluations,completed_populations,completed_pairs,scan_epoch,refinement_epoch,hypercube,mutation_scale,best_fitness,controller_fnv1a32,event\n");
      append(timing_,"completed_evaluations,process_clock_ticks,clocks_per_second\n");
    } catch(...) { if(timing_) fclose(timing_); fclose(trace_); throw; }
  }
  ~Observer() { if(trace_) fclose(trace_); if(timing_) fclose(timing_); }
  unsigned long seed() const { return config_.seed; }
  unsigned long completedPairs() const { return pairs_; }
  bool stopped() const { return !stop_.empty(); }
  bool reserve(unsigned long group) {
    require(!finished_ && reserved_==0,"RESERVATION_STATE");
    if(stopped()) return false;
    require(group==2 || group==40,"GROUP_SIZE");
    if(group>config_.limit-completed_) { stop_="reserve_stop"; return false; }
    reserved_=group; return true;
  }
  void returnedEvaluation() {
    require(!finished_ && !stopped() && reserved_>0 && completed_<config_.limit,"UNRESERVED_EVALUATION");
    ++completed_; --reserved_;
  }
  void boundary(bool population,const std::vector<double>& values,double fitness,
                long scan,long epoch,double hc,double mutation) {
    require(!finished_ && reserved_==0 && completed_>last_boundary_,"BOUNDARY_STATE");
    require(completed_-last_boundary_==(population?40UL:2UL),"BOUNDARY_COUNT");
    require(finite(fitness) && finite(hc) && finite(mutation),"NONFINITE_STATE");
    if(population) ++populations_; else ++pairs_;
    last_boundary_=completed_;
    const std::string eac=controller(values),hash=digest(eac);
    std::ostringstream s; s.imbue(std::locale::classic()); s<<std::setprecision(17);
    s<<completed_<<','<<populations_<<','<<pairs_<<','<<scan<<','<<epoch<<','<<hc<<','<<mutation<<','<<fitness<<','<<hash<<','<<(population?"population":"pair")<<'\n';
    append(trace_,s.str());
    std::ostringstream timing; timing<<(std::clock()-start_)<<','<<CLOCKS_PER_SEC;
    std::ostringstream row; row<<completed_<<','<<timing.str()<<'\n'; append(timing_,row.str());
    while(next_<config_.checkpoints.size() && config_.checkpoints[next_]<=completed_) {
      unsigned long requested=config_.checkpoints[next_++]; actual_.push_back(completed_);
      if(config_.capture) {
        std::ostringstream base; base<<config_.output<<"/checkpoint_"<<std::setw(6)<<std::setfill('0')<<requested;
        write(base.str()+".eac",eac);
        write(base.str()+".json",metadata(requested,"count_checkpoint",fitness,scan,epoch,hc,mutation,hash));
      }
    }
    if(completed_==config_.limit) stop_="evaluation_limit";
  }
  void finish(bool aborted,bool natural,const std::vector<double>& values,double fitness,
              long scan,long epoch,double hc,double mutation) {
    require(!finished_ && reserved_==0 && completed_==last_boundary_,"FINAL_STATE");
    require(!aborted,"USER_OR_SIMULATION_ABORT");
    if(stop_.empty()) { require(natural,"UNCLASSIFIED_STOP"); stop_="natural_convergence"; }
    require(finite(fitness) && finite(hc) && finite(mutation),"NONFINITE_STATE");
    const std::string eac=controller(values),hash=digest(eac);
    write(config_.output+"/final.eac",eac);
    write(config_.output+"/final.json",metadata(completed_,stop_.c_str(),fitness,scan,epoch,hc,mutation,hash));
    std::ostringstream s;
    s<<"{\n\"schema\":1,\"status\":\"complete\",\"study\":\"count_bank_only\",\"seed\":"<<config_.seed
     <<",\"capture\":"<<(config_.capture?"true":"false")<<",\"evaluation_limit\":"<<config_.limit
     <<",\"completed_evaluations\":"<<completed_<<",\"completed_populations\":"<<populations_
     <<",\"completed_pairs\":"<<pairs_<<",\"stop_reason\":\""<<stop_<<"\",\"controller_fnv1a32\":\""<<hash<<"\",\"checkpoints\":[";
    for(size_t i=0;i<config_.checkpoints.size();++i) {
      if(i) s<<',';
      s<<"{\"requested\":"<<config_.checkpoints[i]<<",\"actual\":";
      if(i<actual_.size()) s<<actual_[i]; else s<<"null";
      s<<",\"captured\":"<<((i<actual_.size() && config_.capture)?"true":"false")<<'}';
    }
    s<<"],\"timing_semantics\":\"process_clock_ticks_platform_dependent_not_campaign_wall_time\"}\n";
    write(config_.output+"/summary.json",s.str()); finished_=true;
  }
};
}
#endif
