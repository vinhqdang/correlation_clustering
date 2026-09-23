// Root lower bounds of the KaPoCE branch-and-bound on an instance read from stdin.
#include <iostream>
#include <chrono>
#include <cluster_editing/exact/instance.h>
#include <cluster_editing/exact/lower_bounds.h>
#include <cluster_editing/exact/star_bound.h>
int main() {
  auto inst = load_exact_instance();
  auto t0 = std::chrono::steady_clock::now();
  int p3 = packing_local_search_bound(inst, INF);
  auto t1 = std::chrono::steady_clock::now();
  int st = star_bound(inst, INF);
  auto t2 = std::chrono::steady_clock::now();
  std::cout << "p3 " << p3 << " " << std::chrono::duration<double>(t1 - t0).count() << "\n";
  std::cout << "star " << st << " " << std::chrono::duration<double>(t2 - t1).count() << "\n";
}
