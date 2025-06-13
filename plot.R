#!/usr/bin/env Rscript
# 输出一系列图片和用于合成图片的数据文件
# plot.R <raw_bin file> <platform>
# v0.5e


library(tools)
library(stringi)
library(ggrepel)
args = commandArgs(trailingOnly = TRUE)
raw_bin_file = args[1]
platform = args[2]

#raw_bin_file = 'salus.queryname_95.rawbin.txt'
#platform = 'Saluscall'


output_file = sprintf('%s.RDS', platform)

raw.df = read_tsv(raw_bin_file, comment = '#', show_col_types = F) %>% dplyr::filter(base_called != 'N', base_real != 'N', cycle <= 150)


qual_class <- function(qual) {
   if (qual <= 11) {return('Low')}

   if (12 <= qual & qual <= 25) {return('Medium')}

   if (26 <= qual) {return('High')}
}
raw.df$qual_class = sapply(raw.df$qual, FUN = qual_class)

result.list = list()

# ================================================== summary table ==================================================
total.df = raw.df %>%
   group_by(base_real) %>%
   summarise(bases.mismatch = sum(count[base_called != base_real]), bases.total = sum(count))

temp.df = raw.df %>%
   group_by(base_real, qual_class) %>%
   summarise(mismatch = sum(count[base_called != base_real]))


summary.df = merge(temp.df, total.df, by = 'base_real', all.x = T)
summary.df %<>% mutate(ratio = round(mismatch / bases.mismatch * 100, 2),  error.rate = signif(mismatch / bases.total * 1000, 3))
summary.df %<>% mutate(base_real = factor(base_real, levels = c('A', 'T', 'C', 'G')),
                       qual_class = factor(qual_class, levels = c('Low', 'Medium', 'High')),
                       platform = platform)
summary.df %<>% arrange(base_real, qual_class)

result.list$summary = summary.df
print(summary.df)
write_tsv(summary.df, file = sprintf('summary.%s.tsv', platform))

temp.float = signif(unique(summary.df$bases.mismatch / summary.df$bases.total * 1000), 4)
message('Error rate for A/T/C/G (‰): ')
print(temp.float)

temp.float = signif(sum(summary.df$bases.mismatch) / sum(summary.df$bases.total) * 1000, 4)
message('Average Error rate (‰): ', temp.float)
result.list$error_rate.float = temp.float


# ================================================== reportQ 比例图 ==================================================

temp.df = raw.df %>% group_by(qual) %>% summarise(total = sum(count))
temp.int = sum(temp.df$total)
temp.df %<>% mutate(percent = total / temp.int * 100)

plot <- ggplot(temp.df) +
   geom_bar(aes(x = qual, y = percent, group = qual, fill = qual), position = "stack", stat = "identity", width = 1) +
   geom_text(aes(x = qual, y = percent + 1, label = round(percent, 2)), size = 2) +
   theme(legend.position = 'none',
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         axis.title.x = element_text(size = 10),
         axis.text.x = element_text(size = 10)) +
   xlab('Quality') +
   ylab('Percent of total bases (%)') +
   scale_x_continuous(limits = c(0, 60), expand = expansion(add = c(0, 0)), breaks = seq(0, 60, 2)) +
   scale_y_continuous(limits = c(0, 100), expand = expansion(add = c(0, 0)), breaks = seq(0, 100, 5)) +
   ggtitle(sprintf('Quality percent distribution (%s)', platform))
#print(plot)

ggsave(sprintf('quality.%s.png', platform), plot)
result.list$quality.plot = plot
result.list$quality = temp.df



# ================================================== reportQ ~ readQ plot ==================================================
POINT_SIZE_SHIFT.int = 10

data.df = raw.df %>%
   group_by(qual) %>%
   summarise(mismatch = sum(count[base_called != base_real]),total = sum(count))

data.df %<>% dplyr::filter(total >= 1000, mismatch >= 10) %>%
   mutate(realQ = -log10(mismatch / total) * 10, platform = platform) %>%
   mutate(percentage = total / sum(data.df$total)) %>%
   #mutate(point.size = (log2(percentage) + POINT_SIZE_SHIFT.int))
   mutate(point.size = percentage ^ 0.5 * 20)


temp.df = data.df %>% arrange(desc(percentage)) %>% head(3) %>% mutate(label = sprintf('%s%%', signif(percentage * 100, 3))) %>% select(qual, label)
data.df = merge(data.df, temp.df, by = 'qual', all.x = T)
data.df$label[is.na(data.df$label)] <- ''
print(data.df)

plot <- ggplot() +
   geom_point(data.df, mapping = aes(x = qual, y = realQ, color = point.size, size = point.size), alpha = 0.5) +
   geom_smooth(data.df, mapping = aes(x = qual, y = realQ), span = 0.3, se = F) +
   geom_abline(slope = 1, intercept = 0 , linetype = "dashed", color = "red", alpha = 0.75) +
   scale_x_continuous(limits = c(6, NA), expand = expansion(add = c(0, 0)), breaks = seq(6, 100, 2)) +
   scale_y_continuous(limits = c(6, NA), expand = expansion(add = c(0, 0)), breaks = seq(6, 100, 2)) +
   theme(legend.position = 'none', plot.title = element_text(hjust = 0.5)) +
   ggtitle(sprintf('Q~Q plot (%s)', platform)) +
   xlab('Platform output Q score') +
   ylab('Real Q score') +
   scale_size_identity() +
   ggrepel::geom_text_repel(data = data.df, mapping = aes(x = qual, y = realQ, label = label), nudge_x = 2, nudge_y = 4)
#print(plot)
result.list$qq.ungroup.df = data.df
result.list$qq.loess.upgroup.plot = plot
ggsave(sprintf('qq.loess.ungroup.plot.%s.png', platform), plot)


# 按上游碱基质量分组
data.df = raw.df %>%
   group_by(qual_upstream, qual) %>%
   summarise(mismatch = sum(count[base_called != base_real]),total = sum(count))

data.df %<>% dplyr::filter(total >= 1000, mismatch >= 10, qual_upstream != -1) %>%
   mutate(realQ = -log10(mismatch / total) * 10, platform = platform) %>%
   mutate(percentage = total / sum(data.df$total)) %>%
   #mutate(point.size = (log2(percentage) + POINT_SIZE_SHIFT.int))
   mutate(point.size = percentage ^ 0.5 * 20)

temp.df = data.df %>% arrange(desc(percentage)) %>% head(3) %>% mutate(label = sprintf('%s%%', signif(percentage * 100, 3))) %>% select(qual_upstream, qual, label)
data.df = merge(data.df, temp.df, by = c('qual_upstream', 'qual'), all.x = T)
data.df$label[is.na(data.df$label)] <- ''
print(data.df)

#
plot <- ggplot() +
   geom_point(data.df, mapping = aes(x = qual, y = realQ, color = point.size, size = point.size), alpha = 0.5) +
   geom_smooth(data.df, mapping = aes(x = qual, y = realQ, group = qual_upstream, color = qual_upstream), method = 'lm', se = F) +
   geom_abline(slope = 1, intercept = 0 , linetype = "dashed", color = "red", alpha = 0.75) +
   scale_x_continuous(limits = c(6, 60), expand = expansion(add = c(0, 0)), breaks = seq(6, 60, 2)) +
   scale_y_continuous(limits = c(6, 60), expand = expansion(add = c(0, 0)), breaks = seq(6, 60, 2)) +
   theme(legend.position = 'none', plot.title = element_text(hjust = 0.5)) +
   ggtitle(sprintf('Q~Q plot (%s)', platform)) +
   xlab('Platform output Q score') +
   ylab('Real Q score') +
   scale_size_identity() +
   ggrepel::geom_text_repel(data = data.df, mapping = aes(x = qual, y = realQ, label = label), nudge_x = 2, nudge_y = 4)
#print(plot)
ggsave(sprintf('qq.lm.plot.%s.png', platform), plot)
result.list$qq.df = data.df
result.list$qq.lm.plot = plot

#
plot <- ggplot() +
   geom_point(data.df, mapping = aes(x = qual, y = realQ, color = point.size, size = point.size), alpha = 0.5) +
   geom_smooth(data.df, mapping = aes(x = qual, y = realQ, group = qual_upstream, color = qual_upstream), span = 0.3, se = F) +
   geom_abline(slope = 1, intercept = 0 , linetype = "dashed", color = "red", alpha = 0.75) +
   scale_x_continuous(limits = c(6, 60), expand = expansion(add = c(0, 0)), breaks = seq(6, 60, 2)) +
   scale_y_continuous(limits = c(6, 60), expand = expansion(add = c(0, 0)), breaks = seq(6, 60, 2)) +
   theme(legend.position = 'none', plot.title = element_text(hjust = 0.5)) +
   ggtitle(sprintf('Q~Q plot (%s)', platform)) +
   xlab('Platform output Q score') +
   ylab('Real Q score') +
   scale_size_identity() +
   ggrepel::geom_text_repel(data = data.df, mapping = aes(x = qual, y = realQ, label = label), nudge_x = 2, nudge_y = 4)
#print(plot)
ggsave(sprintf('qq.loess.plot.%s.png', platform), plot)
result.list$qq.loess.plot = plot



# 按上游碱基种类质量分组
data.df = raw.df %>%
   dplyr::filter(!is.na(base_upstream)) %>%
   group_by(base_upstream, base_called, qual) %>%
   summarise(mismatch = sum(count[base_called != base_real]),total = sum(count))

data.df %<>% dplyr::filter(total >= 1000, mismatch >= 10) %>%
   mutate(realQ = -log10(mismatch / total) * 10, platform = platform) %>%
   mutate(percentage = total / sum(data.df$total)) %>%
   #mutate(point.size = (log2(percentage) + POINT_SIZE_SHIFT.int))
   mutate(point.size = percentage ^ 0.5 * 20)

temp.df = data.df %>% arrange(desc(percentage)) %>% head(3) %>% mutate(label = sprintf('%s%%', signif(percentage * 100, 3))) %>% select(base_upstream, base_called, qual, label)
data.df = merge(data.df, temp.df, by = c('base_upstream', 'base_called', 'qual'), all.x = T)
data.df$label[is.na(data.df$label)] <- ''
data.df %<>% mutate(dimer = paste0(base_upstream, base_called)) %>% select(dimer, qual, realQ, label, point.size, platform)
print(data.df)
result.list$qq.base.df = data.df




# ================================================== cycle ~ realQ / mismatch rate plot==================================================
# 按预估Q值分组
temp.df = raw.df %>%
   group_by(cycle, qual) %>%
   summarise(mismatch = sum(count[base_called != base_real]), total = sum(count))


temp.df %<>% mutate(plot.mismatch = mismatch  + 1, plot.total = total + 2)
temp.df %<>% mutate(realQ = -log10(plot.mismatch / plot.total) * 10, platform = platform)
temp.fa = factor(temp.df$qual, levels = sort(unique(temp.df$qual)))
temp.df$qual = temp.fa
result.list$cycle.q.df = temp.df

temp.tb = table(temp.df$qual)
if (all(!is.na(temp.tb[c('11', '25', '37')])) & all(temp.tb[c('11', '25', '37')] >= 75)) {
   AVAILABLE.QUAL = c(11, 25, 37)
} else {
   qual.df = temp.df %>% group_by(qual) %>% summarise(total.bases = sum(total)) %>% arrange(desc(total.bases))  # 取bases最多的三个qual
   AVAILABLE.QUAL = qual.df %>% pull(qual) %>% head(3)
   AVAILABLE.QUAL = as.numeric(as.vector(AVAILABLE.QUAL))
}

temp.df %<>% dplyr::filter(qual %in% AVAILABLE.QUAL, total >= 100)
plot <- ggplot() +
   geom_point(temp.df, mapping = aes(x = cycle, y = realQ, group = qual, color = qual), alpha = 0.5, size = 3) +
   geom_smooth(temp.df, mapping = aes(x = cycle, y = realQ, group = qual, color = qual),  span = 0.8, se = F) +
   theme_classic() +
   theme(legend.position = 'bottom',
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         axis.title.x = element_text(size = 10),
         axis.text.x = element_text(size = 10)) +
   xlab('Cycle') +
   ylab('Real Q score') +
   scale_x_continuous(limits = c(0, NA), breaks = seq(0, 150, 10)) +
   scale_y_continuous(limits = c(0, NA), breaks = seq(0, 100, 5)) +
   ggtitle(sprintf('Q~cycle plot (%s)', platform)) +
   geom_hline(yintercept = AVAILABLE.QUAL, linetype = "dashed", color = "red")
#scale_colour_manual(values = c('11' = 'orangered4', '25' = 'orange2', '37' = 'lightgreen'))
#print(plot)
ggsave(sprintf('cycle.q.plot.%s.png', platform), plot)
result.list$cycle.q.plot = plot


# 按 real base分组
total.df = raw.df %>%
   group_by(base_real, cycle) %>%
   summarise(total = sum(count))

mismatch.df = raw.df %>%
   dplyr::filter(base_called != base_real) %>%
   group_by(base_real, cycle) %>%
   summarise(mismatch.count = sum(count))

data.df = merge(mismatch.df, total.df, by = c('base_real', 'cycle'), all.y = T)
data.df %<>% mutate(error.rate = mismatch.count / total * 1000) %>%
   dplyr::filter(!is.na(error.rate))


total.df = raw.df %>%
   group_by(cycle) %>%
   summarise(total = sum(count))

mismatch.df = raw.df %>%
   dplyr::filter(base_called != base_real) %>%
   group_by(cycle) %>%
   summarise(mismatch.count = sum(count))

overall.df = merge(mismatch.df, total.df, by = c('cycle'), all.y = T)
overall.df %<>% mutate(error.rate = mismatch.count / total * 1000, base_real = 'Overall') %>%
   dplyr::filter(!is.na(error.rate))

# temp.df = raw.df %>%
#    dplyr::filter(base_called != base_real) %>%
#    mutate(mismatch = str_c(base_real, '>', base_called)) %>%
#    group_by(mismatch, cycle) %>%
#    summarise(mismatch.count = sum(total))
data.df = rbind(data.df, overall.df)
data.df %<>% mutate(Real_base = factor(base_real, levels = c('A', 'T', 'C', 'G', 'Overall')), platform = platform)

y_uplimit.float = quantile(data.df$error.rate, 0.95) * 2
plot <- ggplot() +
   geom_point(data.df, mapping = aes(x = cycle, y = error.rate, group = Real_base, color = Real_base), alpha = 0.5) +
   geom_smooth(data.df, mapping = aes(x = cycle, y = error.rate, group = Real_base, color = Real_base), span = 0.5, se = F) +
   scale_x_continuous(limits = c(0, max(data.df$cycle) * 1.05), expand = expansion(add = c(0, 0)), breaks = seq(0, 200, 10)) +
   scale_y_continuous(limits = c(0, y_uplimit.float), expand = expansion(add = c(0, 0)), breaks = seq(0, 10, 0.2)) +
   theme_classic() +
   theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5)) +
   ggtitle(sprintf('Cycle ~ Error rate plot (%s)', platform)) +
   xlab('Cycle') +
   ylab('Error rate (‰)') +
   scale_colour_manual(values = c('A' = 'green', 'T' = 'red', 'C' = 'blue', 'G' = 'orange', 'Overall' = 'black'))
#print(plot)
ggsave(sprintf('cycle.error.rate.plot.%s.png', platform), plot)
result.list$cycle.error.rate.plot = plot
result.list$cycle.error.rate.df = data.df

# 不限制y轴比例
plot <- ggplot() +
   geom_point(data.df, mapping = aes(x = cycle, y = error.rate, group = Real_base, color = Real_base), alpha = 0.5) +
   geom_smooth(data.df, mapping = aes(x = cycle, y = error.rate, group = Real_base, color = Real_base), span = 0.5, se = F) +
   scale_x_continuous(limits = c(0, max(data.df$cycle) * 1.05), expand = expansion(add = c(0, 0)), breaks = seq(0, 200, 10)) +
   scale_y_continuous(limits = c(0, NA), expand = expansion(add = c(0, 0)), breaks = seq(0, 10, 0.2)) +
   theme_classic() +
   theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5)) +
   ggtitle(sprintf('Cycle ~ Error rate plot (%s)', platform)) +
   xlab('Cycle') +
   ylab('Error rate (‰)') +
   scale_colour_manual(values = c('A' = 'green', 'T' = 'red', 'C' = 'blue', 'G' = 'orange', 'Overall' = 'black'))
#print(plot)
ggsave(sprintf('cycle.error.rate.plot2.%s.png', platform), plot)
result.list$cycle.error.rate.plot2 = plot




# 按R1 R2 分组
total.df = raw.df %>%
   group_by(read_order, cycle) %>%
   summarise(mismatch = sum(count[base_called != base_real]), total = sum(count))
total.df %<>% mutate(error.rate = mismatch / total * 1000, read_order = ifelse(read_order == 'f', yes = '1st', no = '2rd'))
y_uplimit.float = quantile(total.df$error.rate, 0.95) * 2

temp.df = total.df %>%
   group_by(read_order) %>%
   summarise(mismatch = sum(mismatch), total = sum(total))
temp.df %<>% mutate(error.rate = mismatch / total * 1000)

plot <- ggplot() +
   geom_point(total.df, mapping = aes(x = cycle, y = error.rate, group = read_order, color = read_order), alpha = 0.5) +
   geom_smooth(total.df, mapping = aes(x = cycle, y = error.rate, group = read_order, color = read_order), span = 0.3, se = F) +
   scale_x_continuous(limits = c(0, max(total.df$cycle) * 1.05), expand = expansion(add = c(0, 0)), breaks = seq(0, 200, 10)) +
   scale_y_continuous(limits = c(0, y_uplimit.float), expand = expansion(add = c(0, 0)), breaks = seq(0, 10, 0.2)) +
   geom_hline(yintercept = temp.df$error.rate[1], color = 'cornflowerblue', linetype = 'dashed') +
   geom_hline(yintercept = temp.df$error.rate[2], color = 'lightpink', linetype = 'dashed') +
   annotate(geom = 'text', x = c(7, 7), y = temp.df$error.rate + 0.02, label = signif(temp.df$error.rate, 3), size = 3) +
   theme_classic() +
   theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5)) +
   ggtitle(sprintf('Cycle ~ Error rate plot (%s)', platform)) +
   xlab('Cycle') +
   ylab('Error rate (‰)') +
   scale_colour_manual(values = c('1st' = 'cornflowerblue', '2rd' = 'lightpink'))
#print(plot)
result.list$cycle.error.rate.read_order.df = total.df
ggsave(sprintf('cycle.error.rate.plot3.%s.png', platform), plot)
result.list$cycle.error.rate.plot3 = plot



# 按 real base (high quality mismatch only)分组
# total.df = raw.df %>%
#    dplyr::filter(qual_class == 'High') %>%
#    group_by(base_real, cycle) %>%
#    summarise(total = sum(count))
#
# mismatch.df = raw.df %>%
#    dplyr::filter(base_called != base_real, qual_class == 'High') %>%
#    group_by(base_real, cycle) %>%
#    summarise(mismatch.count = sum(count))
#
# data.df = merge(mismatch.df, total.df, by = c('base_real', 'cycle'), all.y = T)
# data.df %<>% mutate(error.rate = mismatch.count / total * 1000) %>%
#    dplyr::filter(!is.na(error.rate))
#
#
# total.df = raw.df %>%
#    dplyr::filter(qual_class == 'High') %>%
#    group_by(cycle) %>%
#    summarise(total = sum(count))
#
# mismatch.df = raw.df %>%
#    dplyr::filter(base_called != base_real, , qual_class == 'High') %>%
#    group_by(cycle) %>%
#    summarise(mismatch.count = sum(count))
#
# overall.df = merge(mismatch.df, total.df, by = c('cycle'), all.y = T)
# overall.df %<>% mutate(error.rate = mismatch.count / total * 1000, base_real = 'Overall') %>%
#    dplyr::filter(!is.na(error.rate))
#
# # temp.df = raw.df %>%
# #    dplyr::filter(base_called != base_real) %>%
# #    mutate(mismatch = str_c(base_real, '>', base_called)) %>%
# #    group_by(mismatch, cycle) %>%
# #    summarise(mismatch.count = sum(total))
# data.df = rbind(data.df, overall.df)
# data.df %<>% mutate(Real_base = factor(base_real, levels = c('A', 'T', 'C', 'G', 'Overall')), platform = platform)
# y_uplimit.float = quantile(data.df$error.rate, 0.95) * 2
#
# plot <- ggplot() +
#    geom_point(data.df, mapping = aes(x = cycle, y = error.rate, group = Real_base, color = Real_base), alpha = 0.5) +
#    geom_smooth(data.df, mapping = aes(x = cycle, y = error.rate, group = Real_base, color = Real_base), span = 0.85, se = F) +
#    scale_x_continuous(limits = c(0, max(data.df$cycle) * 1.05), expand = expansion(add = c(0, 0)), breaks = seq(0, 200, 10)) +
#    scale_y_continuous(limits = c(0, y_uplimit.float), expand = expansion(add = c(0, 0)), breaks = seq(0, 10, 0.01)) +
#    theme_classic() +
#    theme(legend.position = 'bottom', plot.title = element_text(hjust = 0.5)) +
#    ggtitle(sprintf('Cycle ~ Error rate plot (High quality only) (%s)', platform)) +
#    xlab('Cycle') +
#    ylab('Error rate (‰)') +
#    scale_colour_manual(values = c('A' = 'green', 'T' = 'red', 'C' = 'blue', 'G' = 'orange', 'Overall' = 'black'))
# #print(plot)
#
# #ggsave(sprintf('cycle.error.rate.high.plot.%s.png', platform), plot)
# #result.list$cycle.error.rate.high.plot = plot
# #result.list$cycle.error.rate.high.df = data.df

# ================================================== cycle ~ real base/mismatch ratio图 ==================================================
# 纵坐标为4种real base在每个cycle内部的占比
total.df = raw.df %>%
   dplyr::filter(base_real != base_called) %>%
   group_by(cycle, base_real) %>%
   summarise(total = sum(count))

temp.df = total.df %>%
   group_by(cycle) %>%
   summarise(total.in.cycle = sum(total))

total.df = merge(total.df, temp.df, by = 'cycle', all.x = T)
total.df %<>% mutate(ratio = total / total.in.cycle * 100)
plot <- ggplot(total.df) +
   geom_bar(aes(x = cycle, y = ratio, group = base_real, fill = base_real), position = "stack", stat = "identity", width = 1) +
   theme_classic() +
   theme(legend.position = 'bottom',
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         axis.title.x = element_text(size = 10),
         axis.text.x = element_text(size = 10)) +
   xlab('Cycle') +
   ylab('Ratio of four real base in each cycle (%)') +
   ggtitle(sprintf('four real base ratio cycle (%s)', platform))
#print(plot)

ggsave(sprintf('real.base.cycle.ratio.plot.%s.png', platform), plot)
result.list$real.base.cycle.ratio.plot = plot
result.list$real.base.cycle.ratio.df = total.df


# 纵坐标为各类mismatch在每个cycle内部的占比
total.df = raw.df %>%
   dplyr::filter(base_real != base_called) %>%
   group_by(cycle, base_real, base_called) %>%
   summarise(total = sum(count))

temp.df = total.df %>%
   group_by(cycle) %>%
   summarise(cycle.total = sum(total))

total.df = merge(total.df, temp.df, by = 'cycle', all.x = T)
total.df %<>% mutate(mismatch = str_c(base_real, '>', base_called),
                     ratio = total / cycle.total * 100)

plot <- ggplot(total.df) +
   #geom_point(total.df, mapping = aes(x = cycle, y = ratio, group = mismatch, color = mismatch), alpha = 0.25, size = 0.5) +
   #geom_smooth(total.df, mapping = aes(x = cycle, y = ratio, group = mismatch, color = mismatch), span = 0.15, se = F) +
   geom_bar(aes(x = cycle, y = ratio, group = mismatch, fill = mismatch), position = "stack", stat = "identity", width = 1) +
   theme_classic() +
   theme(legend.position = 'bottom',
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         axis.title.x = element_text(size = 10),
         axis.text.x = element_text(size = 10)) +
   xlab('Cycle') +
   ylab('Ratio of mismatches in each cycle (%)') +
   #scale_y_continuous(limits = c(0, 100)) +
   ggtitle(sprintf('Mismatch ratio cycle (%s)', platform))
#print(plot)

ggsave(sprintf('mismatch.cycle.ratio.plot.%s.png', platform), plot)
result.list$mismatch.cycle.ratio.plot = plot
result.list$mismatch.cycle.ratio.df = total.df




# 纵坐标为4种real base占全部mismatch的占比
total.df = raw.df %>%
   dplyr::filter(base_real != base_called) %>%
   group_by(cycle, base_real) %>%
   summarise(total = sum(count))
total.mismatch.int = sum(total.df$total)
total.df %<>% mutate(percent = total / total.mismatch.int * 100)

plot <- ggplot() +
   geom_point(total.df, mapping = aes(x = cycle, y = percent, group = base_real, color = base_real), alpha = 0.25, size = 0.5) +
   geom_smooth(total.df, mapping = aes(x = cycle, y = percent, group = base_real, color = base_real), span = 0.15, se = F) +
   theme_classic() +
   theme(legend.position = 'bottom',
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         axis.title.x = element_text(size = 10),
         axis.text.x = element_text(size = 10)) +
   xlab('Cycle') +
   ylab('Percentage of all mismatches (%)') +
   #scale_y_continuous(limits = c(0, 2)) +
   ggtitle(sprintf('Four real base cycle distrbution (%s)', platform))
#print(plot)

ggsave(sprintf('real.base.cycle.distribution.plot.%s.png', platform), plot)
result.list$real.base.cycle.distribution.plot = plot
result.list$real.base.cycle.distribution.df = total.df



# 纵坐标为各类mismatch 占全部mismatch的占比
total.df = raw.df %>%
   dplyr::filter(base_real != base_called) %>%
   group_by(cycle, base_real, base_called) %>%
   summarise(total = sum(count))
total.mismatch.int = sum(total.df$total)
total.df %<>% mutate(mismatch = str_c(base_real, '>', base_called),
                     percent = total / total.mismatch.int * 100)

plot <- ggplot() +
   geom_point(total.df, mapping = aes(x = cycle, y = percent, group = mismatch, color = mismatch), alpha = 0.25, size = 0.5) +
   geom_smooth(total.df, mapping = aes(x = cycle, y = percent, group = mismatch, color = mismatch), span = 0.15, se = F) +
   theme_classic() +
   theme(legend.position = 'bottom',
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         axis.title.x = element_text(size = 10),
         axis.text.x = element_text(size = 10)) +
   xlab('Cycle') +
   ylab('Percentage of all mismatches (%)') +
   #scale_y_continuous(limits = c(0, 2)) +
   ggtitle(sprintf('Mismatch cycle distrbution (%s)', platform))
#print(plot)

ggsave(sprintf('mismatch.cycle.distribution.plot.%s.png', platform), plot)
result.list$mismatch.cycle.distribution.plot = plot
result.list$mismatch.cycle.distribution.df = total.df



# ================================================== mismatch碱基质量分布图 ==================================================
# 横坐标为质量分数Q, 纵坐标为按照4种read base分组的mismatch总数的占比
mismatch.df = raw.df %>%
   dplyr::filter(base_real != base_called) %>%
   group_by(base_real, qual) %>%
   summarise(total = sum(count))
total.int = sum(mismatch.df$total)

mismatch.df %<>% mutate(percent = total / total.int * 100)
plot <- ggplot() +
   geom_point(mismatch.df, mapping = aes(x = qual, y = percent, group = base_real, color = base_real), alpha = 0.25, size = 5) +
   geom_smooth(mismatch.df, mapping = aes(x = qual, y = percent, group = base_real, color = base_real), span = 0.15, se = F) +
   #geom_vline(xintercept = c(12, 26), linetype = "dashed", color = "red") +
   theme_classic() +
   theme(legend.position = 'bottom',
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         axis.title.x = element_text(size = 10),
         axis.text.x = element_text(size = 10)) +
   xlab('Quality') +
   ylab('Percentage of total mismatches (%)') +
   scale_x_continuous(limits = c(0, max(mismatch.df$qual) + 2), expand = expansion(add = c(0, 0)), breaks = seq(0, 100, 2)) +
   scale_y_continuous(limits = c(0, NA), breaks = seq(0, 50, 1)) +
   ggtitle(sprintf('Four real bases Quality distrbution (%s)', platform))
#print(plot)
ggsave(sprintf('real.base.quality.distribution.plot.%s.png', platform), plot)
result.list$real.base.quality.distribution.plot = plot
result.list$real.base.quality.distribution.df = mismatch.df




# 横坐标为质量分数Q, 纵坐标为按照12种mismatch类型分组的mismatch总数的占比
total.df = raw.df %>%
   dplyr::filter(base_real != base_called) %>%
   mutate(mismatch = str_c(base_real, '>', base_called)) %>%
   group_by(mismatch) %>%
   summarise(total = sum(count))

total.int = sum(total.df$total)

mismatch.df = raw.df %>%
   dplyr::filter(base_real != base_called) %>%
   mutate(mismatch = str_c(base_real, '>', base_called))
mismatch.df %<>% group_by(mismatch, qual) %>% summarise(total = sum(count))

data.df = merge(mismatch.df, total.df, by = 'mismatch', all.x = T) %>% mutate(total.z = total.int)
data.df %<>% mutate(ratio.1 = total.x / total.y * 100, ratio.2 = total.x / total.z * 100, platform = platform)

# 纵坐标为各类mismatch内部的占比 (obsoleted)
# plot <- ggplot() +
#    geom_point(data.df, mapping = aes(x = qual, y = ratio.1, group = mismatch, color = mismatch), alpha = 0.25, size = 0.5) +
#    geom_smooth(data.df, mapping = aes(x = qual, y = ratio.1, group = mismatch, color = mismatch),  span = 0.15, se = F) +
#    theme_classic() +
#    theme(legend.position = 'bottom',
#          plot.title = element_text(size = 20, hjust = 0.5, vjust = -1),
#          axis.title.y = element_text(size = 20),
#          axis.text.y = element_text(size = 15),
#          axis.title.x = element_text(size = 20),
#          axis.text.x = element_text(size = 15)) +
#    xlab('Quality') +
#    ylab('Ratio(%)') +
#    scale_y_continuous(limits = c(0, NA)) +
#    ggtitle('Mismatch Bases Quality distrbution')
#print(plot)

plot <- ggplot() +
   geom_point(data.df, mapping = aes(x = qual, y = ratio.2, group = mismatch, color = mismatch), alpha = 0.25, size = 0.5) +
   geom_smooth(data.df, mapping = aes(x = qual, y = ratio.2, group = mismatch, color = mismatch), span = 0.15, se = F) +
   geom_vline(xintercept = c(12, 26), linetype = "dashed", color = "red") +
   theme_classic() +
   theme(legend.position = 'bottom',
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1),
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         axis.title.x = element_text(size = 10),
         axis.text.x = element_text(size = 10)) +
   scale_x_continuous(limits = c(0, max(data.df$qual) + 2), expand = expansion(add = c(0, 0)), breaks = seq(0, 100, 2)) +
   scale_y_continuous(limits = c(0, NA), breaks = seq(0, 50, 0.25)) +
   xlab('Quality') +
   ylab('Percentage of total mismatches (%)') +
   scale_y_continuous(limits = c(0, 2)) +
   ggtitle(sprintf('Mismatch Bases Quality distrbution (%s)', platform))
#print(plot)
ggsave(sprintf('mismatch.quality.distribution.plot.%s.png', platform), plot)
result.list$mismatch.quality.distribution.plot = plot
result.list$mismatch.quality.distribution.df = data.df




# = ================================================== 单碱基偏好 前缀：mismatch.type.base: 错误率柱状图, 比例柱状图 和 grid plot ==================================================
cosine_similarity <- function(A, B) {
   # 计算cosine similarity
   # A和B是长度相等的numeric vector

   if (!is.numeric(A)) {
      message('A must be a numeric vector.')
      return(NULL)
   }

   if (!is.numeric(B)) {
      message('B must be a numeric vector.')
      return(NULL)
   }

   if (length(A) != length(B)) {
      message('A and B must be equal in length.')
      return(NULL)
   }

   numerator = sum(A*B)
   denominator1 = sum(A^2)^0.5
   denominator2 = sum(B^2)^0.5

   return(numerator / (denominator1 * denominator2))
}

total.df = raw.df %>%
   group_by(base_real) %>%
   summarise(total = sum(count))

# 错误率柱状图
mismatch.df = raw.df %>%
   dplyr::filter(base_real != base_called) %>%
   group_by(base_real) %>%
   summarise(mismatch = sum(count))

data.df = merge(mismatch.df, total.df) %>%
   mutate(mismatch.rate = mismatch / total * 1000, platform = platform)
cs.num = cosine_similarity(data.df$mismatch.rate, rep(temp.float, nrow(data.df)))

plot.base.error.rate = ggplot(data = data.df) +
   geom_bar(aes(x = base_real, y = mismatch.rate, fill = base_real), position = position_dodge(0.6), stat = "identity", width = 0.5) +
   geom_hline(yintercept = temp.float, linetype = "dashed", color = "red") +
   scale_y_continuous(limits = c(0, max(data.df$mismatch.rate) * 1.05), expand = expansion(add = c(0, 0))) +
   xlab('Real base') +
   ylab('Error rate (‰)') +
   ggtitle(sprintf('Error rate for each base type (%s)\nCosine similarity:%s', platform, round(cs.num, 3))) +
   theme(legend.position = 'none',
         axis.text.x = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1)) +
   ggplot2::annotate(geom = 'text', x = 0.8, y = temp.float + 0.01, label = sprintf('average: %s', temp.float), size = 3)
#print(plot.base.error.rate)
ggsave(sprintf('mismatch.type.base.error.rate.plot.%s.png', platform), plot.base.error.rate, width = 15, height = 9, units = 'cm')
result.list$mismatch.type.base.error.rate.plot = plot.base.error.rate
result.list$mismatch.type.base.error.rate.df = data.df


# 比例柱状图 (只有base_real)
temp.df = data.df %>%
   select(-platform) %>%
   mutate(mismatch.ratio = mismatch / sum(data.df$mismatch),
                             total.ratio = total / sum(data.df$total)) %>%
   pivot_longer(cols = c(5, 6), names_to = 'type', values_to = 'ratio') %>% select(base_real, type, ratio)
temp.df %<>% mutate(alpha = rep(c(1, 0.95), nrow(data.df)), platform = platform)
cs.num = cosine_similarity(data.df$mismatch / sum(data.df$mismatch), data.df$total / sum(data.df$total))

plot.base.ratio = ggplot(data = temp.df) +
   geom_bar(aes(x = base_real, y = ratio * 100, group = type, fill = base_real, alpha = alpha), position = position_dodge(0.6), stat = "identity", width = 0.5) +
   scale_y_continuous(limits = c(0, max(temp.df$ratio) * 105), expand = expansion(add = c(0, 0))) +
   xlab('Real base') +
   ylab('Ratio (%)') +
   ggtitle(sprintf('Ratio for each base type error on total error (%s)\nCosine similarity:%s', platform, round(cs.num, 3))) +
   theme(legend.position = 'none',
         axis.text.x = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1))
#print(plot.base.ratio)
ggsave(sprintf('mismatch.type.base.ratio.plot.%s.png', platform), plot.base.ratio, width = 13, height = 8, units = 'cm')
result.list$mismatch.type.base.ratio.plot = plot.base.ratio
result.list$mismatch.type.base.ratio.df = temp.df


# 比例柱状图 (附带错误碱基)
mismatch.with.called.df = raw.df %>%
   dplyr::filter(base_real != base_called) %>%
   group_by(base_real, base_called) %>%
   summarise(total = sum(count))

if (nrow(mismatch.with.called.df) != 12) {
   message(platform, '单联碱基偏好 警告: 错误碱基与called base种类总和不为12.')
   print(mismatch.with.called.df, n = 12)
}

temp.df = mismatch.with.called.df %>% mutate(ratio = total / sum(mismatch.with.called.df$total) * 100,
                             label = str_c(base_real, '>', base_called),
                             platform = platform
                             )
plot.base.ratio.with.error = ggplot(data = temp.df) +
   geom_bar(aes(x = label, y = ratio, fill = base_real), position = position_dodge(0.6), stat = "identity", width = 0.5) +
   geom_hline(yintercept = round(1 / 12 * 100, 3), linetype = "dashed", color = "red") +
   scale_y_continuous(limits = c(0, max(temp.df$ratio) * 1.05), expand = expansion(add = c(0, 0))) +
   xlab('Base error types') +
   ylab('Ratio (%)') +
   ggtitle(sprintf('Ratio for each base error type on total error (%s)', platform)) +
   theme(legend.position = 'none',
         axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust = 0.5),
         axis.text.y = element_text(size = 10),
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1))
#print(plot.base.ratio.with.error)
ggsave(sprintf('mismatch.type.base.ratio.with.error.plot.%s.png', platform), plot.base.ratio.with.error, width = 16, height = 8, units = 'cm')
result.list$mismatch.type.base.ratio.with.error.plot = plot.base.ratio.with.error
result.list$mismatch.type.base.ratio.with.error.df = temp.df


# grid plot
data.df = merge(mismatch.with.called.df, mismatch.df, by = 'base_real', all = T)
data.df %<>% mutate(ratio = total / mismatch * 100, platform = platform)
cs.num = round(cosine_similarity(data.df$ratio, rep(1/3 * 100, nrow(data.df))), 3)

plot.base.ratio.grid = ggplot(data = data.df) +
   geom_bar(mapping = aes(x = base_called, y = ratio), position = "dodge", stat = "identity", width = 0.2) +
   facet_wrap(~ base_real, nrow = 2) +
   geom_hline(yintercept = 1/3 * 100, linetype = "dashed", color = "red", linewidth = 0.4) +
   ggtitle(sprintf('Base mismatch type bias (%s)\nCosine similarity:%s', platform, cs.num)) +
   xlab('Called base') +
   ylab('Ratio(%)') +
   theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
#print(plot.base.ratio.grid)
ggsave(sprintf('mismatch.type.base.ratio.grid.plot.%s.png', platform), plot.base.ratio.grid, width = 16, height = 16, units = 'cm')
result.list$mismatch.type.base.ratio.grid.plot = plot.base.ratio.grid
result.list$mismatch.type.base.ratio.grid.df = data.df


# = ================================================== 双联碱基偏好 前缀：mismatch.type.dimer: 错误率柱状图, 比例柱状图 和 grid plot ==================================================
BASE = c('A', 'T', 'C', 'G')
cosine_similarity <- function(A, B) {
   # 计算cosine similarity
   # A和B是长度相等的numeric vector

   if (!is.numeric(A)) {
      message('A must be a numeric vector.')
      return(NULL)
   }

   if (!is.numeric(B)) {
      message('B must be a numeric vector.')
      return(NULL)
   }

   if (length(A) != length(B)) {
      message('A and B must be equal in length.')
      return(NULL)
   }

   numerator = sum(A*B)
   denominator1 = sum(A^2)^0.5
   denominator2 = sum(B^2)^0.5

   return(numerator / (denominator1 * denominator2))
}

total.df = raw.df %>%
   dplyr::filter(!is.na(base_upstream), base_upstream %in% BASE) %>%
   mutate(dimer = str_c(base_upstream, base_real)) %>%
   group_by(dimer) %>%
   summarise(total = sum(count))

if (nrow(total.df) != 16) {
   print(total.df)
   message(platform, ' 双联碱基偏好 警告: 总体二联碱基种类不为16.')
}

mismatch.df = raw.df %>%
   dplyr::filter(!is.na(base_upstream), base_upstream %in% BASE, base_real != base_called) %>%
   mutate(dimer = str_c(base_upstream, base_real)) %>%
   group_by(dimer) %>%
   summarise(mismatch = sum(count))

if (nrow(mismatch.df) != 16) {
   print(mismatch.df)
   message(platform, ' 双联碱基偏好 警告: 错误二联碱基种类不为16.')
}


# 错误率柱状图
data.df = merge(mismatch.df, total.df, by = 'dimer', all = T)
data.df %<>% mutate(mismatch.ratio = mismatch / sum(data.df$mismatch),
                    total.ratio = total / sum(data.df$total),
                    mismatch.rate = mismatch / total * 1000,
                    base = str_sub(dimer, 2, 2),
                    dimer = factor(dimer, levels = stri_reverse(dimer)),
                    platform = platform
                    )

message('Average Error rate (‰): ', temp.float) # 平均错误率（千分之）
cs.num = cosine_similarity(data.df$mismatch.rate, rep(temp.float, nrow(data.df)))
plot.dimer.error.rate = ggplot(data = data.df) +
   geom_bar(aes(x = dimer, y = mismatch.rate, fill = base), position = position_dodge(0.6), stat = "identity", width = 0.5) +
   geom_hline(yintercept = temp.float, linetype = "dashed", color = "red") +
   scale_y_continuous(limits = c(0, max(data.df$mismatch.rate) * 1.05), expand = expansion(add = c(0, 0))) +
   xlab('Dimer types') +
   ylab('Error rate (‰)') +
   ggtitle(sprintf('Error rate for each dimer type (%s)\n Cosine similarity:%s', platform, round(cs.num, 3))) +
   theme(legend.position = 'none',
         axis.text.x = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1)) +
   ggplot2::annotate(geom = 'text', x = 2, y = temp.float + 0.01, label = sprintf('average: %s', temp.float), size = 3)
#print(plot.dimer.error.rate)
ggsave(sprintf('mismatch.type.dimer.error.rate.plot.%s.png', platform), plot.dimer.error.rate, width = 15, height = 6, units = 'cm')
result.list$mismatch.type.dimer.error.rate.plot = plot.dimer.error.rate
result.list$mismatch.type.dimer.error.rate.df = data.df

# 比例柱状图 (只二联碱基)
temp.df = data.df %>% pivot_longer(cols = c(4,5), names_to = 'type', values_to = 'ratio') %>% select(dimer, type, ratio)
temp.df %<>% mutate(alpha = rep(c(1, 0.95), nrow(data.df)),
                    base = str_sub(dimer, 2, 2),
                    dimer = factor(dimer, levels = unique(stri_reverse(dimer))),
                    platform = platform
                    )
cs.num = cosine_similarity(data.df$mismatch.ratio, data.df$total.ratio)

plot.dimer.ratio = ggplot(data = temp.df) +
   geom_bar(aes(x = dimer, y = ratio * 100, group = type, fill = base, alpha = alpha), position = position_dodge(0.6), stat = "identity", width = 0.5) +
   scale_y_continuous(limits = c(0, max(temp.df$ratio) * 105), expand = expansion(add = c(0, 0))) +
   xlab('Dimer types') +
   ylab('Ratio (%)') +
   ggtitle(sprintf('Ratio for each dimer type error on total error (%s)\n Cosine similarity:%s', platform, round(cs.num, 3))) +
   theme(legend.position = 'none',
         axis.text.x = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1))
#print(plot.dimer.ratio)
ggsave(sprintf('mismatch.type.dimer.ratio.plot.%s.png', platform), plot.dimer.ratio, width = 30, height = 8, units = 'cm')
result.list$mismatch.type.dimer.ratio.plot = plot.dimer.ratio
result.list$mismatch.type.dimer.ratio.df = temp.df


# 比例柱状图 (附带错误碱基)
mismatch.with.called.df = raw.df %>%
   dplyr::filter(!is.na(base_upstream), base_upstream %in% BASE, base_real != base_called) %>%
   mutate(dimer = str_c(base_upstream, base_real)) %>%
   group_by(dimer, base_called) %>%
   summarise(total = sum(count))

data.df = merge(mismatch.with.called.df, mismatch.df, by = 'dimer', all = T)
if (nrow(data.df) != 48) {
   print(data.df)
   message(platform, ' 双联碱基偏好 警告: 错误二联碱基与called base种类总和不为48.')
}

temp.df = data.df %>% mutate(ratio = total / sum(data.df$total) * 100,
                             label = str_c(dimer, '>', base_called))
seq.df = data.frame(str_extract_all(temp.df$label, '[ATCG]', simplify = T))
seq.df %<>% arrange(X2, X1, X3)
seq.c = str_c(seq.df$X1, seq.df$X2, '>', seq.df$X3)
temp.df %<>% mutate(label = factor(label, levels = seq.c),
                    base = str_sub(dimer, 2, 2),
                    platform = platform
                    )

plot.dimer.ratio.with.error = ggplot(data = temp.df) +
   geom_bar(aes(x = label, y = ratio, fill = base), position = position_dodge(0.6), stat = "identity", width = 0.5) +
   geom_hline(yintercept = round(1 / 48 * 100, 3), linetype = "dashed", color = "red") +
   scale_y_continuous(limits = c(0, max(temp.df$ratio) * 1.05), expand = expansion(add = c(0, 0))) +
   xlab('Dimer error types') +
   ylab('Ratio (%)') +
   ggtitle(sprintf('Ratio for each dimer error type on total error (%s)', platform)) +
   theme(legend.position = 'none',
         axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust = 0.5),
         axis.text.y = element_text(size = 10),
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1))
#print(plot.dimer.ratio.with.error)
ggsave(sprintf('mismatch.type.dimer.ratio.with.error.plot.%s.png', platform), plot.dimer.ratio.with.error, width = 30, height = 8, units = 'cm')
result.list$mismatch.type.dimer.ratio.with.error.plot = plot.dimer.ratio.with.error
result.list$mismatch.type.dimer.ratio.with.error.df = temp.df


# grid plot
data.df %<>% mutate(ratio = total / mismatch * 100, platform = platform)
cs.num = round(cosine_similarity(data.df$ratio, rep(1/3 * 100, nrow(data.df))), 3)

plot.dimer.ratio.grid = ggplot(data = data.df) +
   geom_bar(mapping = aes(x = base_called, y = ratio), position = "dodge", stat = "identity", width = 0.2) +
   facet_wrap(~ dimer, nrow = 4) +
   geom_hline(yintercept = 1/3 * 100, linetype = "dashed", color = "red", linewidth = 0.4) +
   ggtitle(sprintf('Dimer mismatch type bias (%s)\nCosine similarity:%s', platform, cs.num)) +
   xlab('Called base for the last base in dimer') +
   ylab('Ratio(%)') +
   theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
#print(plot.dimer.ratio.grid)
ggsave(sprintf('mismatch.type.dimer.ratio.grid.plot.%s.png', platform), plot.dimer.ratio.grid, width = 16, height = 20, units = 'cm')
result.list$mismatch.type.dimer.ratio.grid.plot = plot.dimer.ratio.grid
result.list$mismatch.type.dimer.ratio.grid.df = data.df

# = ================================================== 三联碱基偏好 前缀：mismatch.type.trimer: 错误率柱状图, 比例柱状图 和 grid plot==================================================
BASE = c('A', 'T', 'C', 'G')
cosine_similarity <- function(A, B) {
   # 计算cosine similarity
   # A和B是长度相等的numeric vector

   if (!is.numeric(A)) {
      message('A must be a numeric vector.')
      return(NULL)
   }

   if (!is.numeric(B)) {
      message('B must be a numeric vector.')
      return(NULL)
   }

   if (length(A) != length(B)) {
      message('A and B must be equal in length.')
      return(NULL)
   }

   numerator = sum(A*B)
   denominator1 = sum(A^2)^0.5
   denominator2 = sum(B^2)^0.5

   return(numerator / (denominator1 * denominator2))
}

total.df = raw.df %>%
   dplyr::filter(!is.na(base_upstream_2), !is.na(base_upstream), base_upstream_2 %in% BASE, base_upstream %in% BASE) %>%
   mutate(trimer = str_c(base_upstream_2, base_upstream, base_real)) %>%
   group_by(trimer) %>%
   summarise(total = sum(count))

if (nrow(total.df) != 64) {
   print(total.df)
   message(platform, ' 三联碱基偏好 警告: 总体三联碱基种类不为64.')
}

mismatch.df = raw.df %>%
   dplyr::filter(!is.na(base_upstream_2), !is.na(base_upstream), base_upstream_2 %in% BASE, base_upstream %in% BASE, base_real != base_called) %>%
   mutate(trimer = str_c(base_upstream_2, base_upstream, base_real)) %>%
   group_by(trimer) %>%
   summarise(mismatch = sum(count))

if (nrow(mismatch.df) != 64) {
   print(mismatch.df)
   message(platform, ' 三联碱基偏好 警告: 错误三联碱基种类不为64.')
}


# 错误率柱状图
data.df = merge(mismatch.df, total.df, by = 'trimer', all = T)
data.df %<>% mutate(mismatch.ratio = mismatch / sum(data.df$mismatch),
                    total.ratio = total / sum(data.df$total),
                    mismatch.rate = mismatch / total * 1000,
                    base = str_sub(trimer, 3, 3),
                    trimer = factor(trimer, levels = stri_reverse(trimer)),
                    platform = platform
)

message('Average Error rate (‰): ', temp.float) # 平均错误率（千分之）
cs.num = cosine_similarity(data.df$mismatch.rate, rep(temp.float, nrow(data.df)))
plot.trimer.error.rate = ggplot(data = data.df) +
   geom_bar(aes(x = trimer, y = mismatch.rate, fill = base), position = position_dodge(0.6), stat = "identity", width = 0.5) +
   geom_hline(yintercept = temp.float, linetype = "dashed", color = "red") +
   scale_y_continuous(limits = c(0, max(data.df$mismatch.rate) * 1.05), expand = expansion(add = c(0, 0))) +
   xlab('Trimer types') +
   ylab('Error rate (‰)') +
   ggtitle(sprintf('Error rate for each trimer type (%s)\nCosine similarity:%s', platform, round(cs.num, 3))) +
   theme(legend.position = 'none',
         axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust = 0.5),
         axis.text.y = element_text(size = 10),
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1)) +
   ggplot2::annotate(geom = 'text', x = 6, y = temp.float + 0.01, label = sprintf('average: %s', temp.float), size = 3)
#print(plot.trimer.error.rate)
ggsave(sprintf('mismatch.type.trimer.error.rate.plot.%s.png', platform), plot.trimer.error.rate, width = 30, height = 8, units = 'cm')
result.list$mismatch.type.trimer.error.rate.plot = plot.trimer.error.rate
result.list$mismatch.type.trimer.error.rate.df = data.df

# 比例柱状图 (只三联碱基)
temp.df = data.df %>% pivot_longer(cols = c(4,5), names_to = 'type', values_to = 'ratio') %>% select(trimer, type, ratio)
temp.df %<>% mutate(alpha = rep(c(1, 0.95), nrow(data.df)),
                    base = str_sub(trimer, 3, 3),
                    trimer = factor(trimer, levels = unique(stri_reverse(trimer))),
                    platform = platform
)
cs.num = cosine_similarity(data.df$mismatch.ratio, data.df$total.ratio)

plot.trimer.ratio = ggplot(data = temp.df) +
   geom_bar(aes(x = trimer, y = ratio * 100, group = type, fill = base, alpha = alpha), position = position_dodge(0.6), stat = "identity", width = 0.5) +
   scale_y_continuous(limits = c(0, max(temp.df$ratio) * 105), expand = expansion(add = c(0, 0))) +
   xlab('Dimer types') +
   ylab('Ratio (%)') +
   ggtitle(sprintf('Ratio for each trimer type error on total error (%s)\nCosine similarity:%s', platform, round(cs.num, 3))) +
   theme(legend.position = 'none',
         axis.text.x = element_text(size = 10, angle = 45, vjust = 0.5, hjust = 0.5),
         axis.text.y = element_text(size = 10),
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1))
#print(plot.trimer.ratio)
ggsave(sprintf('mismatch.type.trimer.ratio.plot.%s.png', platform), plot.trimer.ratio, width = 30, height = 8, units = 'cm')
result.list$mismatch.type.trimer.ratio.plot = plot.trimer.ratio
result.list$mismatch.type.trimer.ratio.df = temp.df

# 比例柱状图 (附带错误碱基)
mismatch.with.called.df = raw.df %>%
   dplyr::filter(!is.na(base_upstream_2), !is.na(base_upstream), base_upstream_2 %in% BASE, base_upstream %in% BASE, base_real != base_called) %>%
   mutate(trimer = str_c(base_upstream_2, base_upstream, base_real)) %>%
   group_by(trimer, base_called) %>%
   summarise(total = sum(count))

data.df = merge(mismatch.with.called.df, mismatch.df, by = 'trimer', all.x = T)
if (nrow(data.df) != 192) {
   print(data.df)
   message(platform, ' 三联碱基偏好 警告: 错误三联碱基与called base种类总和不为192.')
}

temp.df = data.df %>% mutate(ratio = total / sum(data.df$total) * 100,
                             label = str_c(trimer, '>', base_called))
seq.df = data.frame(str_extract_all(temp.df$label, '[ATCG]', simplify = T))
seq.df %<>% arrange(X3, X4, X2, X1)
seq.c = str_c(seq.df$X1, seq.df$X2, seq.df$X3, '>', seq.df$X4)
temp.df %<>% mutate(label = factor(label, levels = seq.c),
                    base = str_sub(trimer, 3, 3),
                    platform = platform
)

plot.trimer.ratio.with.error = ggplot(data = temp.df) +
   geom_bar(aes(x = label, y = ratio, fill = base), position = position_dodge(0.6), stat = "identity", width = 0.5) +
   geom_hline(yintercept = round(1 / 192 * 100, 3), linetype = "dashed", color = "red") +
   scale_y_continuous(limits = c(0, max(temp.df$ratio) * 1.05), expand = expansion(add = c(0, 0))) +
   xlab('Trimer error types') +
   ylab('Ratio (%)') +
   ggtitle(sprintf('Ratio for each trimer error type on total error (%s)', platform)) +
   theme(legend.position = 'none',
         axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 0.5),
         axis.text.y = element_text(size = 10),
         plot.title = element_text(size = 10, hjust = 0.5, vjust = -1))
#print(plot.trimer.ratio.with.error)
ggsave(sprintf('mismatch.type.trimer.ratio.with.error.plot.%s.png', platform), plot.trimer.ratio.with.error, width = 60, height = 8, units = 'cm')
result.list$mismatch.type.trimer.ratio.with.error.plot = plot.trimer.ratio.with.error
result.list$mismatch.type.trimer.ratio.with.error.df = temp.df


# grid plot
data.df %<>% mutate(ratio = total / mismatch * 100,
                    trimer = factor(trimer, levels = unique(stri_reverse(trimer))),
                    platform = platform
                    )
cs.num = round(cosine_similarity(data.df$ratio, rep(1/3 * 100, nrow(data.df))), 3)

plot.trimer.ratio.grid = ggplot(data = data.df) +
   geom_bar(mapping = aes(x = base_called, y = ratio), position = "dodge", stat = "identity", width = 0.2) +
   facet_wrap(~ trimer, nrow = 8) +
   geom_hline(yintercept = 1/3 * 100, linetype = "dashed", color = "red", linewidth = 0.4) +
   ggtitle(sprintf('Trimer mismatch type bias (%s)\nCosine similarity:%s', platform, cs.num)) +
   xlab('Called base for the last base in trimer') +
   ylab('Ratio(%)') +
   theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))
#print(plot.trimer.ratio.grid)
ggsave(sprintf('mismatch.type.trimer.ratio.grid.plot.%s.png', platform), plot.trimer.ratio.grid, width = 16, height = 20, units = 'cm')
result.list$mismatch.type.trimer.ratio.grid.plot = plot.trimer.ratio.grid
result.list$mismatch.type.trimer.ratio.grid.df = data.df


# = ==================================================  每种 mismatch的 prefix base 占比图  ==================================================

data.df = raw.df %>% dplyr::filter(base_real != base_called, !is.na(base_upstream), base_upstream != 'N') %>%
   mutate(mismatch = str_c(base_real, '->', base_called)) %>%
   group_by(base_upstream, mismatch) %>%
   summarise(total = sum(count))

temp.df = data.df %>%
   group_by(mismatch) %>%
   summarise(mismatch.total = sum(total))


data.df = merge(data.df, temp.df, by = 'mismatch', all.x = T)
data.df %<>% mutate(percentage = total / mismatch.total * 100)



plot.prefix.base.percent = ggplot(data = data.df) +
   geom_bar(mapping = aes(x = base_upstream, y = percentage), position = position_dodge(0.6), stat = "identity", width = 0.4) +
   facet_wrap(~ mismatch, ncol = 3) +
   geom_hline(yintercept = 1/4 * 100, linetype = "dashed", color = "red", linewidth = 0.4) +
   xlab('Prefix base') +
   ylab('Percentage(%)') +
   theme(legend.position = 'bottom',
         plot.title = element_text(hjust = 0.5),
         axis.title.x = element_text(size = 6),
         axis.title.y = element_text(size = 6),
         strip.text.x = element_text(size = 6)
   )
#print(plot.prefix.base.percent)

ggsave(sprintf('mismatch.type.prefix.base.plot.%s.png', platform), plot.prefix.base.percent, width = 15, height = 15, units = 'cm')
result.list$mismatch.type.prefix.base.plot = plot.prefix.base.percent
result.list$mismatch.type.prefix.base.df = data.df




saveRDS(result.list, file = output_file)
message(output_file, ' Done.')
