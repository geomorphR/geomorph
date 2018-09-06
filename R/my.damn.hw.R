
unicef = read.csv("https://raw.githubusercontent.com/mlcollyer/chatham.bio532/master/Data/unicef.csv")

bw <- na.omit(unicef$lowbwt)

# Question 2.17

# a) construct a box plot for the percentages of low birthweight patients

boxplot(bw, ylab = "Birth Weight")

# b) Do the data appeared to be skewed?  Yes, they are slightly positively skewed, but this 
# is mostly because of outliers

# c) Yes, two values are larger than 1.5 X IQR, on the positive end
