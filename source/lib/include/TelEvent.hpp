#pragma once

#include <vector>
#include <memory>
#include <cstddef>
#include <algorithm>
namespace altel{

  struct TelMeasRaw{
    union {
      uint64_t index;
      uint16_t loc[4];
      uint16_t uvdcUS[4];
      int16_t  uvdcS[4];
    } data{0};

    TelMeasRaw(uint64_t h)
      :data{ .index = h }{};
    TelMeasRaw(uint16_t u, uint16_t v, uint16_t detN, uint16_t clk)
      :data{ .uvdcUS = {u, v, detN, clk} }{};

    inline bool operator==(const TelMeasRaw &rh) const{
      return data.index == rh.data.index;
    }

    inline bool operator<(const TelMeasRaw &rh) const{
      return data.index < rh.data.index;
    }

    inline const uint64_t& index() const  {return data.index;}
    inline const uint16_t& u() const  {return data.loc[0];}
    inline const uint16_t& v() const  {return data.loc[1];}
    inline const uint16_t& detN() const  {return data.loc[2];}
    inline const uint16_t& clk() const  {return data.loc[3];}

    inline uint64_t& index(){return data.index;}
    inline uint16_t& u() {return data.loc[0];}
    inline uint16_t& v() {return data.loc[1];}
    inline uint16_t& detN() {return data.loc[2];}
    inline uint16_t& clk() {return data.loc[3];}
  };

  struct TelMeasHit{
    uint16_t DN{0};    // detector id
    double PLs[2]{0.0, 0.0};  // local pos
    std::vector<TelMeasRaw> MRs;// measures

    inline const double& u() const  {return PLs[0];}
    inline const double& v() const  {return PLs[1];}
    inline const uint16_t& detN() const  {return DN;}

    inline double& u() {return PLs[0];}
    inline double& v() {return PLs[1];}
    inline uint16_t& detN() {return DN;}

    static std::vector<std::shared_ptr<TelMeasHit>>
    clustering_UVDCus(const std::vector<TelMeasRaw>& measRaws){
      std::vector<std::shared_ptr<TelMeasHit>> measHits;

      auto hit_col_remain= measRaws;
      while(!hit_col_remain.empty()){
        std::vector<TelMeasRaw> hit_col_this_cluster;
        std::vector<TelMeasRaw> hit_col_this_cluster_edge;

        // get first edge seed hit
        // from un-identifed hit to edge hit
        hit_col_this_cluster_edge.push_back(hit_col_remain[0]);
        hit_col_remain.erase(hit_col_remain.begin());

        while(!hit_col_this_cluster_edge.empty()){
          auto ph_e = hit_col_this_cluster_edge[0];
          uint64_t e = ph_e.index();
          uint64_t c = 0x00000001;
          uint64_t r = 0x00010000;
          //  8 sorround hits search
          std::vector<TelMeasRaw> sorround_col
            {e-c+r, e+r, e+c+r,
             e-c  ,      e+c,
             e-c-r, e-r, e+c-r
            };
          for(auto& sr: sorround_col){
            // only search on un-identifed hits
            auto sr_found_it = std::find(hit_col_remain.begin(), hit_col_remain.end(), sr);
            if(sr_found_it != hit_col_remain.end()){
              // move the found sorround hit
              // from un-identifed hit to an edge hit
              hit_col_this_cluster_edge.push_back(sr);
              hit_col_remain.erase(sr_found_it);
            }
          }
          // after sorround search
          // move from edge hit to cluster hit
          hit_col_this_cluster.emplace_back(e);
          hit_col_this_cluster_edge.erase(hit_col_this_cluster_edge.begin());
        }

        //no edge pixel any more, full cluster
        std::shared_ptr<TelMeasHit> measHit(new TelMeasHit);
        measHit->MRs=std::move(hit_col_this_cluster);
        measHits.push_back(measHit);
      }

      for(auto &measHit : measHits){
        uint16_t detN = measHit->MRs[0].detN();
        measHit->DN=detN;

        uint64_t sumRawU=0;
        uint64_t sumRawV=0;
        uint64_t rawN=measHit->MRs.size();
        for(auto &measRaw : measHit->MRs){
          sumRawU+= measRaw.u();
          sumRawV+= measRaw.v();
        }
        constexpr double offsetU = -0.02924 * (1024/2 - 0.5);
        constexpr double offsetV = -0.02688 * (512/2 - 0.5);
        measHit->PLs[0] = double(sumRawU)/rawN*0.02924 + offsetU;
        measHit->PLs[1] = double(sumRawV)/rawN*0.02688 + offsetV;
      }
      return measHits;
    }
  };

  struct TelFitHit{
    uint16_t DN{0};   // detector id
    double PLs[2]{0.0, 0.0}; // local pos
    double DLs[3]{0.0, 0.0, 0.0}; // local dir

    double PGs[3]{0.0, 0.0, 0.0}; // global pos
    double DGs[3]{0.0, 0.0, 0.0}; // global dir

    std::shared_ptr<TelMeasHit> OM; // origin measure hit
  };

  struct TelTrajHit{
    uint16_t DN{0}; // detector id
    std::shared_ptr<TelFitHit>  FH;  // fitted hit
    std::shared_ptr<TelMeasHit> MM; // measure hit, matched

    bool hasFitHit() const{
      if(FH){
        return true;
      }
      else{
        return false;
      }
    }

    bool hasMatchedMeasHit() const{
      if(FH){
        return true;
      }
      else{
        return false;
      }
    }

    bool hasOriginMeasHit() const{
      if( FH && FH->OM){
        return true;
      }
      else{
        return false;
      }
    }
  };

  struct TelTrajectory{
    uint64_t TN{0};                               // trajectories id
    std::vector<std::shared_ptr<TelTrajHit>> THs; // hit
    std::shared_ptr<TelTrajHit> trajHit(size_t id){
      for(const auto& th: THs){
        if (th && th->DN == id){
          return th;
        }
      }
      return nullptr;
    }

    size_t numTrajHit()const{
      return THs.size();
    }

    size_t numFitHit()const{
      size_t n = 0;
      for(const auto& th: THs){
        if ( th && th->FH){
          n++;
        }
      }
      return n;
    }

    size_t numOriginMeasHit()const{
      size_t n = 0;
      for(const auto& th: THs){
        if ( th && th->FH && th->FH->OM){
          n++;
        }
      }
      return n;
    }

    size_t numMatchedMeasHit()const{
      size_t n = 0;
      for(const auto& th: THs){
        if ( th && th->MM){
          n++;
        }
      }
      return n;
    }
  };

  struct TelEvent{
    uint32_t RN{0}; // run id
    uint32_t EN{0}; // event id
    uint16_t DN{0}; // detector/setup id
    uint64_t CK{0}; // timestamp / trigger id;
    std::vector<TelMeasRaw> MRs;// measures raw
    std::vector<std::shared_ptr<TelMeasHit>> MHs;// measure hit
    std::vector<std::shared_ptr<TelTrajectory>> TJs; // trajectories

    TelEvent()
      :RN(0), EN(0), DN(0), CK(0){};

    TelEvent(uint32_t rN, uint32_t eN, uint16_t detN, uint64_t clk)
      :RN(rN), EN(eN), DN(detN), CK(clk){};

  };
}
