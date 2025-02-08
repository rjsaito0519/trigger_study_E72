#include <iostream>
#include <Math/LorentzVector.h>
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>
#include <Rtypes.h>

int main() {
    // using namespace ROOT::Math;
    
    using LorentzVec = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<Double_t>>;

    LorentzVec p1(0.0, 0.0, 2.0, 2.938);
    LorentzVec p2(0.0, 0.0, 0.0, 0.939);
    ROOT::Math::XYZVector v = p1.Vect();
    std::cout << v.X() << ", " << v.Y() << ", " << v.Z() << std::endl;

    LorentzVec Ptot = p1 + p2;

    // ブーストベクトル
    auto beta = Ptot.BoostToCM();
    std::cout << beta.X() << ", " << beta.Y() << ", " << beta.Z() << std::endl;

    // Boost オブジェクトを作成
    ROOT::Math::Boost boostCM(beta);

    // 重心系へ Boost
    LorentzVec p1_cm = boostCM(p1);
    LorentzVec p2_cm = boostCM(p2);

    // 結果を表示
    std::cout << "Lab frame:\n";
    std::cout << "p1: (" << p1.Px() << ", " << p1.Py() << ", " << p1.Pz() << ", " << p1.E() << ")\n";
    std::cout << "p2: (" << p2.Px() << ", " << p2.Py() << ", " << p2.Pz() << ", " << p2.E() << ")\n";

    std::cout << "CM frame:\n";
    std::cout << "p1_cm: (" << p1_cm.Px() << ", " << p1_cm.Py() << ", " << p1_cm.Pz() << ", " << p1_cm.E() << ")\n";
    std::cout << "p2_cm: (" << p2_cm.Px() << ", " << p2_cm.Py() << ", " << p2_cm.Pz() << ", " << p2_cm.E() << ")\n";

    return 0;
}
