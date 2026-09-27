// Root lower bounds of the KaPoCE branch-and-bound on an instance read from stdin.
//   lbounds            prints "p3 LB SECONDS" and "star LB SECONDS"
//   lbounds STARS      computes the star bound only and writes its stars to the
//                      file STARS (one per line, centre first; 0-based vertices)
#include <iostream>
#include <fstream>
#include <chrono>
#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/lower_bounds.h>
#include <cluster_editing/exact/star_bound.h>
int star_bound_dump(const Instance &inst, int limit, std::ostream &out);
int main(int argc, char **argv) {
  auto inst = load_exact_instance();
  if (argc > 1) {
    std::ofstream out(argv[1]);
    auto t1 = std::chrono::steady_clock::now();
    int st = star_bound_dump(inst, INF, out);
    auto t2 = std::chrono::steady_clock::now();
    out.close();
    std::cout << "star " << st << " " << std::chrono::duration<double>(t2 - t1).count() << std::endl;
    return 0;
  }
  auto t0 = std::chrono::steady_clock::now();
  int p3 = packing_local_search_bound(inst, INF);
  auto t1 = std::chrono::steady_clock::now();
  std::cout << "p3 " << p3 << " " << std::chrono::duration<double>(t1 - t0).count() << std::endl;
  int st = star_bound(inst, INF);
  auto t2 = std::chrono::steady_clock::now();
  std::cout << "star " << st << " " << std::chrono::duration<double>(t2 - t1).count() << "\n";
}
