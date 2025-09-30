// std::vector<std::pair<int,int>>
// CellListParallel::cellsWithParity(Parity par) const
// {
//     std::vector<std::pair<int,int>> out;
//     out.reserve((nx_in_ * ny_in_) / 4);

//     // Global offset measured in *interior cells*
//     const int offx = cx_ * nx_in_;
//     const int offy = cy_ * ny_in_;

//     for (int iy = 1; iy <= ny_in_; ++iy) {
//         for (int ix = 1; ix <= nx_in_; ++ix) {
//             const int gx = offx + (ix - 1);
//             const int gy = offy + (iy - 1);
//             const bool ex = (gx % 2) == 0;
//             const bool ey = (gy % 2) == 0;
//             const Parity cur = ex ? (ey ? Parity::EvenEven : Parity::EvenOdd)
//                                   : (ey ? Parity::OddEven  : Parity::OddOdd);
//             if (cur == par) out.emplace_back(ix, iy);
//         }
//     }
//     return out;
// }
