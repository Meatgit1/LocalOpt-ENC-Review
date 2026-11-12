#ifndef INTRANCE_NATIVE_BRIDGE_PROBE_H_INCLUDED
#define INTRANCE_NATIVE_BRIDGE_PROBE_H_INCLUDED

// Isolated ideal-sail propagation diagnostic, not mission optimisation.
// Unset environment: no I/O, altered controls, extra EOM call or physics change.
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <CSC.h>

namespace NativeBridgeProbe
{
  enum Mode { Disabled, Export, Replay, Schedule };
  enum Profile { ShortAct, MercuryBodyRendezvousV1, MercuryBodyRendezvousPrecisionV2,
                 MercuryBodyRendezvousPrecisionV3 };
  struct Command
  {
    double time,clock,cone;
    Command(double t,double a,double b) : time(t),clock(a),cone(b) {}
  };
  inline bool& failure_flag() { static bool failed=false; return failed; }
  inline void mark_failed() { failure_flag()=true; }
  inline void require(bool condition,const char* message)
  {
    if (!condition) { mark_failed(); throw std::runtime_error(std::string("NATIVE_BRIDGE: ")+message); }
  }
  inline bool finite(double x) { return std::isfinite(x)!=0; }
  inline bool is_mercury_profile(Profile profile)
  { return profile==MercuryBodyRendezvousV1 || profile==MercuryBodyRendezvousPrecisionV2 ||
           profile==MercuryBodyRendezvousPrecisionV3; }
  inline const char* profile_name(Profile profile)
  {
    if(profile==ShortAct)return "act_short_v1";
    if(profile==MercuryBodyRendezvousV1)return "mercury_body_rendezvous_v1";
    if(profile==MercuryBodyRendezvousPrecisionV2)return "mercury_body_rendezvous_precision_v2";
    require(profile==MercuryBodyRendezvousPrecisionV3,"unknown diagnostic profile");
    return "mercury_body_rendezvous_precision_v3";
  }
  inline bool accepts_mercury_precision(Profile profile,double relative_error,double absolute_error)
  {
    // V1/V2 are frozen. V3 adds a separate, explicit precision-study pair;
    // these are integrator tolerances, not mission/validation acceptance limits.
    if(!finite(relative_error) || !finite(absolute_error) || relative_error!=absolute_error)return false;
    if(profile==MercuryBodyRendezvousV1)return relative_error==1e-6;
    if(profile==MercuryBodyRendezvousPrecisionV2)return relative_error==1e-8 || relative_error==1e-10;
    return profile==MercuryBodyRendezvousPrecisionV3 && (relative_error==1e-12 || relative_error==1e-14);
  }
  inline double parse_number(const std::string& field)
  {
    require(!field.empty(),"empty command field");
    std::istringstream input(field); double value;
    require(static_cast<bool>(input>>value),"invalid command number");
    input>>std::ws;
    require(input.eof(),"trailing characters in command field");
    require(finite(value),"non-finite command number"); return value;
  }
  inline std::vector<double> parse_values(const std::string& line,std::size_t expected)
  {
    require(!line.empty(),"blank command record");
    std::vector<double> values; std::string::size_type start=0;
    for (;;)
    {
      require(values.size()<expected,"too many CSV fields");
      const std::string::size_type comma=line.find(',',start);
      values.push_back(parse_number(line.substr(start,comma==std::string::npos?std::string::npos:comma-start)));
      if (comma==std::string::npos) break;
      start=comma+1;
    }
    require(values.size()==expected,"incorrect CSV field count"); return values;
  }
  inline void check_angles(double a,double b)
  {
    require(finite(a) && a>=-3.14159265358979323846 && a<=3.14159265358979323846,
            "clock outside native [-pi,pi] range");
    require(finite(b) && b>=0 && b<=1.57079632679489661923,
            "ideal-sail cone outside native [0,pi/2] range");
  }
  inline Command parse_command_line(const std::string& line)
  {
    const std::vector<double> v=parse_values(line,3); check_angles(v[1],v[2]);
    return Command(v[0],v[1],v[2]);
  }
  struct ScheduleLimits
  {
    double min_duration,max_duration,max_interval,min_radius,max_radius;
    ScheduleLimits() : min_duration(0),max_duration(0),max_interval(0),min_radius(0),max_radius(0) {}
  };
  inline const char* limits_header()
  { return "min_duration_d,max_duration_d,max_interval_d,min_radius_AU,max_radius_AU"; }
  inline ScheduleLimits parse_schedule_limits(const std::string& line,Profile profile=ShortAct)
  {
    require(profile==ShortAct || is_mercury_profile(profile),"unknown diagnostic profile");
    const std::vector<double> v=parse_values(line,5); ScheduleLimits s;
    s.min_duration=v[0];s.max_duration=v[1];s.max_interval=v[2];s.min_radius=v[3];s.max_radius=v[4];
    if(profile==ShortAct)
      require(s.min_duration>0 && s.max_duration>=s.min_duration && s.max_duration<=10,
              "duration contract outside short diagnostic domain");
    else
      require(s.min_duration>=0.5 && s.max_duration>=s.min_duration && s.max_duration<=600,
              "duration contract outside checked Mercury diagnostic domain");
    require(s.max_interval>0 && s.max_interval<=2,"interval contract outside (0,2] days");
    require(s.min_radius>0 && s.max_radius>s.min_radius && s.max_radius<5,
            "radius contract must remain strictly inside (0,5) AU");
    if(is_mercury_profile(profile))
      require(s.min_radius>=0.2 && s.max_radius<=1.5,"Mercury radius contract outside [0.2,1.5] AU");
    return s;
  }
  inline void validate_schedule(const std::vector<Command>& commands,const ScheduleLimits& s,double epoch)
  {
    require(commands.size()>=2 && commands.size()<=1000,"schedule row count outside [2,1000]");
    require(finite(epoch) && commands.front().time==epoch,"schedule initial epoch differs from fixed reference");
    const double duration=commands.back().time-epoch;
    require(finite(duration) && duration>=s.min_duration && duration<=s.max_duration,
            "schedule duration outside declared limits");
    for (std::size_t j=0;j<commands.size();++j)
    {
      require(finite(commands[j].time),"non-finite command epoch");
      check_angles(commands[j].clock,commands[j].cone);
      if (j) { const double h=commands[j].time-commands[j-1].time;
        require(h>0 && h<=s.max_interval,"schedule intervals must be positive and within declared limit"); }
    }
  }
  inline std::string read_environment(const char* name)
  {
#if defined(_MSC_VER)
    char* allocated=0;std::size_t length=0;
    const int status=::_dupenv_s(&allocated,&length,name);
    std::unique_ptr<char,void(*)(void*)> value(allocated,&std::free);
    require(status==0,"cannot read bridge environment");
    return value.get()==0?std::string():std::string(value.get());
#else
    const char* value=std::getenv(name);
    return value==0?std::string():std::string(value);
#endif
  }
  inline Profile parse_profile(const std::string& value)
  {
    if(value.empty() || value=="act_short_v1")return ShortAct;
    if(value=="mercury_body_rendezvous_v1")return MercuryBodyRendezvousV1;
    if(value=="mercury_body_rendezvous_precision_v2")return MercuryBodyRendezvousPrecisionV2;
    require(value=="mercury_body_rendezvous_precision_v3","unknown INTRANCE_BRIDGE_PROFILE");
    return MercuryBodyRendezvousPrecisionV3;
  }
  struct Probe
  {
    Mode mode;Profile profile;
    bool opened,finished,active,initial_verified,endpoint_recorded,interval_active,reference_set;
    bool target_endpoint_recorded;
    unsigned long updates,rhs_rows,passes_started,passes_completed,completed_intervals;
    unsigned long profile_guard_checks;
    double last_time,interval_start,interval_end,min_seen_radius,max_seen_radius,max_clock_error;
    double reference_mp,reference_srd;
    ScheduleLimits limits;std::vector<Command> imported;std::vector<double> expected_initial;
    std::ofstream commands,states,rhs,initial,intervals,endpoint,target_endpoint;
    Probe() : mode(Disabled),profile(ShortAct),opened(false),finished(false),active(false),initial_verified(false),
      endpoint_recorded(false),interval_active(false),reference_set(false),target_endpoint_recorded(false),updates(0),rhs_rows(0),
      passes_started(0),passes_completed(0),completed_intervals(0),profile_guard_checks(0),last_time(0),interval_start(0),interval_end(0),
      min_seen_radius(1e100),max_seen_radius(0),max_clock_error(0),reference_mp(0),reference_srd(0)
    {
      const std::string value=read_environment("INTRANCE_BRIDGE_MODE"); if (value.empty()) return;
      if (value=="export") mode=Export;else if (value=="replay") mode=Replay;
      else if (value=="schedule") mode=Schedule;
      else require(false,"INTRANCE_BRIDGE_MODE must be export, replay, schedule (or unset)");
      profile=parse_profile(read_environment("INTRANCE_BRIDGE_PROFILE"));
      require(!is_mercury_profile(profile) || mode==Export || mode==Schedule,
              "Mercury profile supports only export or explicit schedule modes");
    }
  private: Probe(const Probe&); Probe& operator=(const Probe&);
  };
  inline Probe& probe() { static Probe value;return value; }
  inline bool configured() { return probe().mode!=Disabled; }
  inline bool enabled() { return configured() && probe().active; }
  inline bool imports_controls() { return probe().mode==Replay || probe().mode==Schedule; }
  inline bool replaying() { return enabled() && imports_controls(); }
  inline bool scheduling() { return enabled() && probe().mode==Schedule; }
  inline bool mercury_profile() { return is_mercury_profile(probe().profile); }
  inline void note_mercury_profile_checked()
  {
    Probe& p=probe();require(enabled() && mercury_profile() && p.profile_guard_checks==p.updates,
                            "Mercury profile guard must run exactly once per command boundary");
    ++p.profile_guard_checks;
  }
  class ReportingPass
  {
    bool engaged,completed;
  public:
    ReportingPass() : engaged(false),completed(false)
    {
      if (!configured()) return;Probe& p=probe();
      require(!failure_flag(),"earlier error before reporting pass");
      require(!p.active && !p.opened && !p.finished && p.passes_started==0,"reporting pass may begin exactly once");
      p.active=true;++p.passes_started;engaged=true;
    }
    void complete()
    {
      if (!engaged) return;Probe& p=probe();
      require(!completed && p.active && !failure_flag(),"invalid reporting pass completion");
      require(p.opened && p.updates>=2 && p.rhs_rows>0,"reporting pass did not export a trajectory");
      require(!imports_controls() || p.updates==p.imported.size(),"reporting pass did not consume full history");
      require(!imports_controls() || p.initial_verified,"initial state was not verified");
      require(p.mode!=Schedule || (p.endpoint_recorded && !p.interval_active && p.completed_intervals+1==p.updates),
              "schedule did not reach and record explicit endpoint");
      require(!mercury_profile() || (p.profile_guard_checks==p.updates &&
              (p.mode!=Schedule || p.target_endpoint_recorded)),"Mercury profile/target endpoint incomplete");
      p.active=false;++p.passes_completed;completed=true;
    }
    ~ReportingPass() { if (engaged && !completed) { mark_failed();probe().active=false; } }
  private: ReportingPass(const ReportingPass&);ReportingPass& operator=(const ReportingPass&);
  };
  inline void open_new(std::ofstream& stream,const char* name)
  {
    std::ifstream old(name);require(!old.good(),"output already exists; use a fresh isolated directory");old.close();
    stream.open(name,std::ios::out);require(stream.good(),"cannot create diagnostic output");
    stream<<std::scientific<<std::setprecision(17);
  }
  inline void strip_cr(std::string& line) { if (!line.empty() && line[line.size()-1]=='\r') line.erase(line.size()-1); }
  inline const char* initial_header()
  { return "MJD_d,r_AU,phi_rad,theta_rad,dr_AU_per_d,dphi_rad_per_d,dtheta_rad_per_d,mp_kg,srd_native,x_AU,y_AU,z_AU,vx_AU_per_d,vy_AU_per_d,vz_AU_per_d,lightness,ideal_sail_flag"; }
  inline void initialise()
  {
    Probe& p=probe();if (p.mode==Disabled || p.opened) return;
    require(p.active && !p.finished,"diagnostic initialisation outside active unfinished reporting pass");
    require(!mercury_profile() || p.profile_guard_checks==1,"Mercury profile not checked before initialisation");
    if (imports_controls())
    {
      std::ifstream input("bridge_input.csv");require(input.good(),"missing bridge_input.csv");std::string line;
      require(static_cast<bool>(std::getline(input,line)),"empty command file");strip_cr(line);
      require(line=="MJD_d,clock_rad,cone_rad","unexpected ideal-sail command header");
      while (std::getline(input,line))
      {
        strip_cr(line);const Command q=parse_command_line(line);
        require(p.imported.empty() || q.time>p.imported.back().time,"command times not strictly increasing");
        require(p.imported.size()<1000,"command file exceeds diagnostic row limit");p.imported.push_back(q);
      }
      require(input.eof() && p.imported.size()>=2,"incomplete command history");
      std::ifstream initial_input("bridge_initial_input.csv");require(initial_input.good(),"missing bridge_initial_input.csv");
      require(static_cast<bool>(std::getline(initial_input,line)),"empty initial-state file");strip_cr(line);
      require(line==initial_header(),"unexpected initial-state header");
      require(static_cast<bool>(std::getline(initial_input,line)),"missing initial-state row");
      p.expected_initial=parse_values(line,17);
      require(!std::getline(initial_input,line) && initial_input.eof(),"extra initial-state records or read error");
      if (p.mode==Schedule)
      {
        std::ifstream limits_input("bridge_schedule_limits.csv");require(limits_input.good(),"missing bridge_schedule_limits.csv");
        require(static_cast<bool>(std::getline(limits_input,line)),"empty schedule limits file");strip_cr(line);
        require(line==limits_header(),"unexpected schedule limits header");
        require(static_cast<bool>(std::getline(limits_input,line)),"missing schedule limits row");
        p.limits=parse_schedule_limits(line,p.profile);
        require(!std::getline(limits_input,line) && limits_input.eof(),"extra schedule limits records or read error");
        validate_schedule(p.imported,p.limits,p.expected_initial[0]);
      }
    }
    open_new(p.commands,"bridge_commands.csv");open_new(p.states,"bridge_states.csv");
    open_new(p.rhs,"bridge_rhs.csv");open_new(p.initial,"bridge_initial.csv");
    p.commands<<"MJD_d,clock_rad,cone_rad\n";
    p.states<<"command_index,MJD_d,r_AU,phi_rad,theta_rad,dr_AU_per_d,dphi_rad_per_d,dtheta_rad_per_d,mp_kg,srd_native,x_AU,y_AU,z_AU,vx_AU_per_d,vy_AU_per_d,vz_AU_per_d,clock_rad,cone_rad,lightness,ideal_sail_flag\n";
    p.rhs<<"rhs_index,command_index,MJD_stage_d,clock_rad,cone_rad,zeta_rad,r_AU,phi_rad,theta_rad,dr_AU_per_d,dphi_rad_per_d,dtheta_rad_per_d,mp_kg,srd_native,yp_r_AU_per_d,yp_phi_rad_per_d,yp_theta_rad_per_d,yp_dr_AU_per_d2,yp_dphi_rad_per_d2,yp_dtheta_rad_per_d2,yp_mp_kg_per_d,yp_srd_native_per_d,a_r_AU_per_d2,a_phi_AU_per_d2,a_theta_AU_per_d2,lightness,ideal_sail_flag\n";
    p.initial<<initial_header()<<'\n';
    if (p.mode==Schedule)
    {
      open_new(p.intervals,"bridge_intervals.csv");open_new(p.endpoint,"bridge_endpoint.csv");
      p.intervals<<"interval_index,start_MJD_d,requested_end_MJD_d,actual_end_MJD_d,duration_d,RK_status,reached_boundary,mp_kg,srd_native,radius_AU\n";
      p.endpoint<<"requested_final_MJD_d,actual_final_MJD_d,duration_d,reached_end,completed_intervals,command_rows,min_sampled_radius_AU,max_sampled_radius_AU,max_boundary_clock_error_d,path_checks_passed,initial_state_verified,x_AU,y_AU,z_AU,vx_AU_per_d,vy_AU_per_d,vz_AU_per_d,final_mp_kg,final_srd_native,diagnostic_only\n";
      if(mercury_profile())
      {
        open_new(p.target_endpoint,"bridge_target_endpoint.csv");
        p.target_endpoint<<"requested_final_MJD_d,actual_final_MJD_d,target_body_id,target_MJD_d,spacecraft_x_AU,spacecraft_y_AU,spacecraft_z_AU,spacecraft_vx_AU_per_d,spacecraft_vy_AU_per_d,spacecraft_vz_AU_per_d,target_x_AU,target_y_AU,target_z_AU,target_vx_AU_per_d,target_vy_AU_per_d,target_vz_AU_per_d,position_residual_AU,velocity_residual_AU_per_d,diagnostic_only\n";
      }
    }
    p.opened=true;
  }
  inline void check_spacecraft(const CSC& sc)
  {
    require(sc.isSolarSail() && sc.isIdealSail(),"only ideal solar sail supported");
    require(sc.getCenterBodyId()==SpaceLib::SUN,"only heliocentric fixture supported");
    require(sc.getPropSys()->getEomCoordinateType()==SpaceLib::POLAR,"only active polar EOM supported");
    require(finite(sc.getLightness()) && sc.getLightness()>0,"invalid ideal-sail lightness");
    check_angles(sc.getClockAngle(),sc.getConeAngle());
  }
  inline void check_state(double epoch,const MathLib::Vector& y)
  {
    if (!enabled()) return;Probe& p=probe();
    require(p.opened && y.size()==8 && finite(epoch),"invalid diagnostic state");
    for (int j=0;j<8;++j) require(finite(y(j)),"non-finite diagnostic state component");
    require(y(0)>0,"nonpositive radius");
    require(std::fabs(std::cos(y(2)))>1e-8,"polar singularity outside diagnostic domain");
    require(std::sqrt(y(4)*y(4)*std::cos(y(2))*std::cos(y(2))+y(5)*y(5))>1e-12,
            "degenerate instantaneous orbit frame");
    if (p.reference_set)
    {
      require(std::fabs(y(6)-p.reference_mp)<=1e-10,"sail auxiliary mass drift");
      require(y(7)>=p.reference_srd-1e-12,"sail dose fell below initial exposure");
    }
    if (y(0)<p.min_seen_radius)p.min_seen_radius=y(0);
    if (y(0)>p.max_seen_radius)p.max_seen_radius=y(0);
    if (scheduling()) require(y(0)>=p.limits.min_radius && y(0)<=p.limits.max_radius,
                             "scheduled state outside declared solar-distance bounds");
  }
  inline void check_schedule_stage(double epoch,const MathLib::Vector& y)
  {
    if (!enabled()) return;Probe& p=probe();
    require(p.updates>0 && p.reference_set,"RHS before initial state/command record");
    if (scheduling()) require(p.interval_active && epoch>=p.interval_start-1e-9 && epoch<=p.interval_end+1e-9,
                             "RHS stage outside active scheduled interval");
    check_state(epoch,y);
  }
  inline bool has_next_interval()
  {
    Probe& p=probe();require(scheduling() && p.opened && p.initial_verified && !p.endpoint_recorded,
                            "schedule queried before initialization or after endpoint");
    return p.updates<p.imported.size();
  }
  inline double begin_interval(const CSC& sc)
  {
    Probe& p=probe();require(has_next_interval() && !p.interval_active && p.updates==p.completed_intervals+1,
                            "invalid schedule interval order");
    check_spacecraft(sc);check_state(sc.MJD(),sc.PolarStatePropellantMassSRDVector());
    require(std::fabs(sc.MJD()-p.last_time)<=1e-9,"spacecraft changed epoch between boundaries");
    p.interval_start=sc.MJD();p.interval_end=p.imported[p.updates].time;
    const double h=p.interval_end-p.interval_start;
    require(finite(h) && h>0 && h<=p.limits.max_interval+1e-9,"invalid next scheduled interval length");
    p.interval_active=true;return h;
  }
  inline void finish_interval(const CSC& sc,int native_status,bool native_done)
  {
    Probe& p=probe();require(scheduling() && p.interval_active,"integration returned without scheduled interval");
    const double error=std::fabs(sc.MJD()-p.interval_end);const bool reached=finite(sc.MJD()) && error<=1e-9;
    const MathLib::Vector y=sc.PolarStatePropellantMassSRDVector();
    p.intervals<<p.completed_intervals<<','<<p.interval_start<<','<<p.interval_end<<','<<sc.MJD()
      <<','<<p.interval_end-p.interval_start<<','<<native_status<<','<<(reached?1:0)<<','<<y(6)<<','<<y(7)<<','<<y(0)<<'\n';
    p.intervals.flush();require(p.intervals.good(),"interval audit write failed");
    require(native_done,"native integrator did not return RK_DONE");require(reached,"scheduled boundary not reached");
    check_spacecraft(sc);check_state(sc.MJD(),y);
    if(error>p.max_clock_error)p.max_clock_error=error;
    ++p.completed_intervals;p.interval_active=false;
  }
  inline void import_command(CSC& sc)
  {
    Probe& p=probe();require(replaying(),"import called outside replay pass");initialise();check_spacecraft(sc);
    require(p.updates<p.imported.size(),"native loop requested too many commands");
    const Command& q=p.imported[p.updates];require(std::fabs(q.time-sc.MJD())<=1e-9,"native epoch differs from command");
    sc.performSteering(q.clock,q.cone);
    require(sc.getClockAngle()==q.clock && sc.getConeAngle()==q.cone,"native setter changed imported command");
  }
  inline void record_command(const CSC& sc)
  {
    Probe& p=probe();if(!enabled())return;initialise();check_spacecraft(sc);
    require(!mercury_profile() || p.profile_guard_checks==p.updates+1,"missing Mercury guard at command boundary");
    const double t=sc.MJD();require(finite(t) && p.updates<1000 && (p.updates==0 || t>p.last_time),
                                   "repeated/unordered command boundary or row limit exceeded");
    const MathLib::Vector y=sc.PolarStatePropellantMassSRDVector();const MathLib::Vector xyz=sc.CState::Cart();
    require(y.size()==8 && xyz.size()==6,"unexpected native state dimensions");check_state(t,y);
    if(p.updates==0)
    {
      std::vector<double> values;values.push_back(t);
      for(int j=0;j<8;++j)values.push_back(y(j));for(int j=0;j<6;++j)values.push_back(xyz(j));
      values.push_back(sc.getLightness());values.push_back(sc.isIdealSail()?1.0:0.0);
      require(values.size()==17,"unexpected initial-state field count");
      for(std::size_t j=0;j<values.size();++j)
      {
        require(finite(values[j]),"non-finite initialized state");
        if(imports_controls()) require(p.expected_initial.size()==values.size() && values[j]==p.expected_initial[j],
                                       "initialized spacecraft/epoch/model differs from exported reference");
        if(j)p.initial<<',';p.initial<<values[j];
      }
      p.initial<<'\n';p.initial.flush();require(p.initial.good(),"initial-state write failed");
      p.initial_verified=imports_controls();p.reference_mp=y(6);p.reference_srd=y(7);p.reference_set=true;
    }
    p.commands<<t<<','<<sc.getClockAngle()<<','<<sc.getConeAngle()<<'\n';p.states<<p.updates<<','<<t;
    for(int j=0;j<8;++j)p.states<<','<<y(j);
    for(int j=0;j<6;++j){require(finite(xyz(j)),"non-finite Cartesian state");p.states<<','<<xyz(j);}
    p.states<<','<<sc.getClockAngle()<<','<<sc.getConeAngle()<<','<<sc.getLightness()<<",1\n";
    p.last_time=t;++p.updates;p.commands.flush();p.states.flush();
    require(p.commands.good() && p.states.good(),"command/state write failed");
  }
  inline void record_sail_rhs(double t,const MathLib::Vector& y,const MathLib::Vector& yp,
                             const CSC& sc,double ar,double aphi,double atheta,double zeta)
  {
    Probe& p=probe();if(!enabled())return;
    require(p.opened && p.updates>0 && !p.finished && p.rhs_rows<250000,"EOM export outside active history or row limit");
    require(y.size()==8 && yp.size()==8,"unexpected EOM state dimensions");check_spacecraft(sc);
    require(finite(t) && finite(ar) && finite(aphi) && finite(atheta) && finite(zeta),"non-finite EOM diagnostic");
    require(yp(6)==0.0,"ideal-sail mass derivative is not zero");
    p.rhs<<p.rhs_rows<<','<<p.updates-1<<','<<t<<','<<sc.getClockAngle()<<','<<sc.getConeAngle()<<','<<zeta;
    for(int j=0;j<8;++j){require(finite(y(j)),"non-finite EOM state");p.rhs<<','<<y(j);}
    for(int j=0;j<8;++j){require(finite(yp(j)),"non-finite native derivative");p.rhs<<','<<yp(j);}
    p.rhs<<','<<ar<<','<<aphi<<','<<atheta<<','<<sc.getLightness()<<",1\n";
    ++p.rhs_rows;if(p.rhs_rows%256==0)p.rhs.flush();require(p.rhs.good(),"EOM export write failed");
  }
  inline void record_endpoint(const CSC& sc)
  {
    Probe& p=probe();require(scheduling() && p.initial_verified && !p.endpoint_recorded && !p.interval_active,
                            "invalid explicit endpoint lifecycle");
    require(p.updates==p.imported.size() && p.completed_intervals+1==p.updates,"endpoint without full schedule");
    require(std::fabs(sc.MJD()-p.imported.back().time)<=1e-9,"actual endpoint differs from requested final epoch");
    check_spacecraft(sc);const MathLib::Vector y=sc.PolarStatePropellantMassSRDVector();check_state(sc.MJD(),y);
    const MathLib::Vector xyz=sc.CState::Cart();require(xyz.size()==6,"invalid endpoint dimension");
    p.endpoint<<p.imported.back().time<<','<<sc.MJD()<<','<<sc.MJD()-p.expected_initial[0]<<",1,"
      <<p.completed_intervals<<','<<p.updates<<','<<p.min_seen_radius<<','<<p.max_seen_radius<<','<<p.max_clock_error<<",1,1";
    for(int j=0;j<6;++j){require(finite(xyz(j)),"non-finite endpoint state");p.endpoint<<','<<xyz(j);}
    p.endpoint<<','<<y(6)<<','<<y(7)<<",1\n";p.endpoint.flush();require(p.endpoint.good(),"endpoint write failed");
    p.endpoint_recorded=true;
    std::cout<<"NATIVE_SCHEDULE_REACHED_END requested_MJD="<<std::setprecision(17)<<p.imported.back().time
      <<" actual_MJD="<<sc.MJD()<<" continuous_from_original_initial_state=yes shooting_resets=0 diagnostic_only=yes"<<std::endl;
  }
  inline void record_target_endpoint(const CSC& sc,int target_body,double target_epoch,
                                     const MathLib::Vector& target)
  {
    Probe& p=probe();require(scheduling() && mercury_profile() && p.endpoint_recorded && !p.target_endpoint_recorded,
                            "invalid Mercury target endpoint lifecycle");
    require(finite(target_epoch) && target_epoch==sc.MJD(),"target epoch differs from actual endpoint");
    const MathLib::Vector xyz=sc.CState::Cart();require(xyz.size()==6 && target.size()==6,"invalid encounter state dimensions");
    double dr2=0,dv2=0;
    p.target_endpoint<<p.imported.back().time<<','<<sc.MJD()<<','<<target_body<<','<<target_epoch;
    for(int j=0;j<6;++j){require(finite(xyz(j)),"non-finite encounter spacecraft state");p.target_endpoint<<','<<xyz(j);}
    for(int j=0;j<6;++j)
    {
      require(finite(target(j)),"non-finite native target state");p.target_endpoint<<','<<target(j);
      const double difference=xyz(j)-target(j);if(j<3)dr2+=difference*difference;else dv2+=difference*difference;
    }
    require(finite(dr2) && finite(dv2),"non-finite encounter residual");
    p.target_endpoint<<','<<std::sqrt(dr2)<<','<<std::sqrt(dv2)<<",1\n";p.target_endpoint.flush();
    require(p.target_endpoint.good(),"target endpoint write failed");p.target_endpoint_recorded=true;
  }
  inline void finish()
  {
    Probe& p=probe();if(p.mode==Disabled)return;require(!failure_flag(),"earlier diagnostic/native error; refusing completion");
    require(!p.active && p.passes_started==1 && p.passes_completed==1,"exactly one completed reporting pass required");
    require(p.opened && p.updates>=2 && p.rhs_rows>0,"diagnostic did not run");
    require(!imports_controls() || (p.updates==p.imported.size() && p.initial_verified),"incomplete imported history/initial comparison");
    require(p.mode!=Schedule || (p.endpoint_recorded && !p.interval_active && p.completed_intervals+1==p.updates),
            "scheduled endpoint not recorded");
    require(!mercury_profile() || (p.profile_guard_checks==p.updates &&
            (p.mode!=Schedule || p.target_endpoint_recorded)),"Mercury profile/target endpoint incomplete at finish");
    p.commands.flush();p.states.flush();p.rhs.flush();p.initial.flush();
    require(p.commands.good() && p.states.good() && p.rhs.good() && p.initial.good(),"final export flush failed");
    p.commands.close();p.states.close();p.rhs.close();p.initial.close();
    if(p.mode==Schedule){p.intervals.flush();p.endpoint.flush();require(p.intervals.good() && p.endpoint.good(),"schedule flush failed");
      p.intervals.close();p.endpoint.close();}
    if(mercury_profile() && p.mode==Schedule){p.target_endpoint.flush();
      require(p.target_endpoint.good(),"target endpoint flush failed");p.target_endpoint.close();}
    p.finished=true;
    std::cout<<"NATIVE_BRIDGE_COMPLETE propulsion=IDEAL_SOLARSAIL mode="
      <<(p.mode==Export?"export":(p.mode==Schedule?"schedule":"replay"))<<" commands="<<p.updates<<" rhs_rows="<<p.rhs_rows
      <<" initial_state_verified="<<(p.initial_verified?"yes":"exported")<<" imported_commands="<<p.imported.size()
      <<" nn_steering_updates="<<(p.mode==Export?p.updates:0)<<" replay_NN_csv_placeholders="<<(imports_controls()?"zero":"none")
      <<" recorded_pass=reporting passes_started="<<p.passes_started<<" passes_completed="<<p.passes_completed
      <<" preparatory_NN_evaluation=unchanged diagnostic_only=yes"<<std::endl;
    if(p.mode==Schedule)std::cout<<"NATIVE_SCHEDULE_COMPLETE reached_end=1 completed_intervals="<<p.completed_intervals
      <<" terminal_command=metadata continuous_initialization_count=1 shooting_resets=0 diagnostic_only=yes"<<std::endl;
    if(mercury_profile())std::cout<<"NATIVE_BRIDGE_PROFILE_COMPLETE profile="<<profile_name(p.profile)<<" profile_guard_checks="
      <<p.profile_guard_checks<<" target_endpoint_recorded="<<(p.target_endpoint_recorded?1:0)<<" diagnostic_only=yes"<<std::endl;
  }
}
#endif
