class Cell
{
 public:
  Cell() = default;
  Cell(int cellId, bool hg, float a, float t, short s) : absId(cellId), isHG(hg), amp(a), time(t), status(s) {}
  ~Cell() = default;
  /// \brief Comparison oparator, based on absId and HG/LG
  /// \param another cell
  /// \return result of comparison: first absId, if eq., first HG, then LG
  inline bool operator<(const Cell& other) const
  {
    if (absId == other.absId) {
      return isHG;
    } else {
      return absId < other.absId;
    }
  }
  /// \brief Comparison oparator, based on absId and HG/LG
  /// \param another cell
  /// \return result of comparison: first absId, if eq., first HG, then LG
  inline bool operator>(const Cell& other) const
  {
    if (absId == other.absId) {
      return !isHG;
    } else {
      return absId > other.absId;
    }
  }

 public:
  short absId = 0;  // abs ID
  bool isHG = 0;    // caloFlag
  short status = 0; // overflow, no Time etc
  float amp = 0;    // Extracted amplitide
  float time = 0;   // extracted time

  ClassDefNV(Cell, 1);
};
