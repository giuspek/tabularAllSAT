(set-logic QF_UF)
(declare-fun A0 () Bool)
(declare-fun A1 () Bool)
(declare-fun A2 () Bool)
(declare-fun A3 () Bool)
(declare-fun A4 () Bool)
(declare-fun A5 () Bool)
(declare-fun A6 () Bool)
(declare-fun A7 () Bool)
(declare-fun A8 () Bool)
(declare-fun A9 () Bool)
(declare-fun A10 () Bool)
(declare-fun A11 () Bool)
(declare-fun A12 () Bool)
(declare-fun A13 () Bool)
(declare-fun A14 () Bool)
(declare-fun A15 () Bool)
(declare-fun A16 () Bool)
(declare-fun A17 () Bool)
(declare-fun A18 () Bool)
(declare-fun A19 () Bool)
(declare-fun A20 () Bool)
(declare-fun A21 () Bool)
(declare-fun A22 () Bool)
(declare-fun A23 () Bool)
(declare-fun A24 () Bool)
(assert (let ((.def_0 (not A18))) (let ((.def_1 (not A15))) (let ((.def_2 (or .def_1 .def_0))) (let ((.def_3 (not .def_2))) (let ((.def_4 (or A23 A9))) (let ((.def_5 (not .def_4))) (let ((.def_6 (or .def_5 .def_3))) (let ((.def_7 (not .def_6))) (let ((.def_8 (not A13))) (let ((.def_9 (or A21 .def_8))) (let ((.def_10 (not .def_9))) (let ((.def_11 (or A2 .def_8))) (let ((.def_12 (not .def_11))) (let ((.def_13 (or .def_12 .def_10))) (let ((.def_14 (and .def_13 .def_7))) (let ((.def_15 (not A16))) (let ((.def_16 (or A5 .def_15))) (let ((.def_17 (not A8))) (let ((.def_18 (or .def_17 .def_15))) (let ((.def_19 (and .def_18 .def_16))) (let ((.def_20 (not A12))) (let ((.def_21 (not A10))) (let ((.def_22 (and .def_21 .def_20))) (let ((.def_23 (not .def_22))) (let ((.def_24 (not A22))) (let ((.def_25 (not A11))) (let ((.def_26 (and .def_25 .def_24))) (let ((.def_27 (and .def_26 .def_23))) (let ((.def_28 (or .def_27 .def_19))) (let ((.def_29 (not .def_28))) (let ((.def_30 (and .def_29 .def_14))) (let ((.def_31 (not A5))) (let ((.def_32 (and A14 .def_31))) (let ((.def_33 (not .def_32))) (let ((.def_34 (not A17))) (let ((.def_35 (not A7))) (let ((.def_36 (or .def_35 .def_34))) (let ((.def_37 (not .def_36))) (let ((.def_38 (= .def_37 .def_33))) (let ((.def_39 (= A15 .def_35))) (let ((.def_40 (not A20))) (let ((.def_41 (or A6 .def_40))) (let ((.def_42 (and .def_41 .def_39))) (let ((.def_43 (and .def_42 .def_38))) (let ((.def_44 (= A9 A21))) (let ((.def_45 (not A0))) (let ((.def_46 (and .def_45 .def_24))) (let ((.def_47 (or .def_46 .def_44))) (let ((.def_48 (or A2 A21))) (let ((.def_49 (not A19))) (let ((.def_50 (or .def_49 A2))) (let ((.def_51 (or .def_50 .def_48))) (let ((.def_52 (not .def_51))) (let ((.def_53 (= .def_52 .def_47))) (let ((.def_54 (not .def_53))) (let ((.def_55 (= .def_54 .def_43))) (let ((.def_56 (and .def_55 .def_30))) (let ((.def_57 (and A18 A10))) (let ((.def_58 (not A21))) (let ((.def_59 (or .def_58 A9))) (let ((.def_60 (not .def_59))) (let ((.def_61 (= .def_60 .def_57))) (let ((.def_62 (not .def_61))) (let ((.def_63 (and .def_0 .def_15))) (let ((.def_64 (or .def_49 .def_8))) (let ((.def_65 (or .def_64 .def_63))) (let ((.def_66 (or .def_65 .def_62))) (let ((.def_67 (= .def_45 .def_15))) (let ((.def_68 (not A14))) (let ((.def_69 (and .def_68 .def_20))) (let ((.def_70 (or .def_69 .def_67))) (let ((.def_71 (or A15 .def_17))) (let ((.def_72 (and A0 .def_8))) (let ((.def_73 (not .def_72))) (let ((.def_74 (and .def_73 .def_71))) (let ((.def_75 (not .def_74))) (let ((.def_76 (or .def_75 .def_70))) (let ((.def_77 (or .def_76 .def_66))) (let ((.def_78 (not .def_77))) (let ((.def_79 (= A0 .def_17))) (let ((.def_80 (not .def_79))) (let ((.def_81 (and .def_49 .def_21))) (let ((.def_82 (not .def_81))) (let ((.def_83 (and .def_82 .def_80))) (let ((.def_84 (and .def_17 .def_25))) (let ((.def_85 (not .def_84))) (let ((.def_86 (or .def_1 A13))) (let ((.def_87 (not .def_86))) (let ((.def_88 (and .def_87 .def_85))) (let ((.def_89 (not .def_88))) (let ((.def_90 (or .def_89 .def_83))) (let ((.def_91 (or A20 A23))) (let ((.def_92 (not .def_91))) (let ((.def_93 (not A23))) (let ((.def_94 (and A6 .def_93))) (let ((.def_95 (not .def_94))) (let ((.def_96 (or .def_95 .def_92))) (let ((.def_97 (or A20 A16))) (let ((.def_98 (and A8 .def_17))) (let ((.def_99 (or .def_98 .def_97))) (let ((.def_100 (not .def_99))) (let ((.def_101 (or .def_100 .def_96))) (let ((.def_102 (or .def_101 .def_90))) (let ((.def_103 (and .def_102 .def_78))) (let ((.def_104 (not .def_103))) (let ((.def_105 (= .def_104 .def_56))) (let ((.def_106 (and .def_0 A2))) (let ((.def_107 (not .def_106))) (let ((.def_108 (or .def_25 A1))) (let ((.def_109 (not .def_108))) (let ((.def_110 (and .def_109 .def_107))) (let ((.def_111 (not .def_110))) (let ((.def_112 (and A17 A7))) (let ((.def_113 (not .def_112))) (let ((.def_114 (and A6 .def_8))) (let ((.def_115 (not .def_114))) (let ((.def_116 (and .def_115 .def_113))) (let ((.def_117 (not .def_116))) (let ((.def_118 (and .def_117 .def_111))) (let ((.def_119 (or .def_21 A8))) (let ((.def_120 (not A6))) (let ((.def_121 (and .def_120 .def_120))) (let ((.def_122 (not .def_121))) (let ((.def_123 (or .def_122 .def_119))) (let ((.def_124 (and A18 .def_25))) (let ((.def_125 (or A16 A12))) (let ((.def_126 (not .def_125))) (let ((.def_127 (or .def_126 .def_124))) (let ((.def_128 (not .def_127))) (let ((.def_129 (or .def_128 .def_123))) (let ((.def_130 (and .def_129 .def_118))) (let ((.def_131 (not .def_130))) (let ((.def_132 (or A9 A15))) (let ((.def_133 (not .def_132))) (let ((.def_134 (not A24))) (let ((.def_135 (not A1))) (let ((.def_136 (and .def_135 .def_134))) (let ((.def_137 (or .def_136 .def_133))) (let ((.def_138 (not .def_137))) (let ((.def_139 (or A9 .def_40))) (let ((.def_140 (not .def_139))) (let ((.def_141 (= .def_8 A2))) (let ((.def_142 (not .def_141))) (let ((.def_143 (or .def_142 .def_140))) (let ((.def_144 (not .def_143))) (let ((.def_145 (= .def_144 .def_138))) (let ((.def_146 (or A16 A2))) (let ((.def_147 (not A2))) (let ((.def_148 (or .def_24 .def_147))) (let ((.def_149 (or .def_148 .def_146))) (let ((.def_150 (not .def_149))) (let ((.def_151 (or .def_134 A6))) (let ((.def_152 (and A3 A15))) (let ((.def_153 (not .def_152))) (let ((.def_154 (and .def_153 .def_151))) (let ((.def_155 (not .def_154))) (let ((.def_156 (and .def_155 .def_150))) (let ((.def_157 (or .def_156 .def_145))) (let ((.def_158 (and .def_157 .def_131))) (let ((.def_159 (not .def_158))) (let ((.def_160 (or .def_135 .def_135))) (let ((.def_161 (and A15 .def_21))) (let ((.def_162 (and .def_161 .def_160))) (let ((.def_163 (not .def_162))) (let ((.def_164 (or A17 A8))) (let ((.def_165 (not .def_164))) (let ((.def_166 (= .def_15 .def_58))) (let ((.def_167 (not .def_166))) (let ((.def_168 (or .def_167 .def_165))) (let ((.def_169 (and .def_168 .def_163))) (let ((.def_170 (not .def_169))) (let ((.def_171 (or .def_34 A6))) (let ((.def_172 (and .def_35 .def_31))) (let ((.def_173 (or .def_172 .def_171))) (let ((.def_174 (and A2 .def_35))) (let ((.def_175 (= A22 .def_0))) (let ((.def_176 (not .def_175))) (let ((.def_177 (or .def_176 .def_174))) (let ((.def_178 (not .def_177))) (let ((.def_179 (or .def_178 .def_173))) (let ((.def_180 (not .def_179))) (let ((.def_181 (and .def_180 .def_170))) (let ((.def_182 (and .def_8 A21))) (let ((.def_183 (not A9))) (let ((.def_184 (not A3))) (let ((.def_185 (or .def_184 .def_183))) (let ((.def_186 (= .def_185 .def_182))) (let ((.def_187 (not .def_186))) (let ((.def_188 (and .def_40 .def_40))) (let ((.def_189 (not .def_188))) (let ((.def_190 (or .def_45 A9))) (let ((.def_191 (or .def_190 .def_189))) (let ((.def_192 (not .def_191))) (let ((.def_193 (and .def_192 .def_187))) (let ((.def_194 (or A13 A1))) (let ((.def_195 (and A16 A4))) (let ((.def_196 (not .def_195))) (let ((.def_197 (and .def_196 .def_194))) (let ((.def_198 (= A16 .def_45))) (let ((.def_199 (not .def_198))) (let ((.def_200 (and A4 A2))) (let ((.def_201 (not .def_200))) (let ((.def_202 (and .def_201 .def_199))) (let ((.def_203 (and .def_202 .def_197))) (let ((.def_204 (and .def_203 .def_193))) (let ((.def_205 (and .def_204 .def_181))) (let ((.def_206 (and .def_205 .def_159))) (let ((.def_207 (not .def_206))) (let ((.def_208 (or .def_207 .def_105))) .def_208))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
(check-sat)
